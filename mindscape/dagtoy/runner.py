# mindscape/dagtoy/runner.py
import importlib
import argparse
import hashlib
from graphviz import Digraph
from pathlib import Path
import yaml

def topological_sort(graph):
    visited, order = set(), []
    def visit(node):
        if node not in visited:
            for dep in graph.get(node, []):
                visit(dep)
            visited.add(node)
            order.append(node)
    for node in graph:
        visit(node)
    return order

def snake_case(name):
    import re
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()

def compute_meta_hash(config_path):
    with open(config_path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()

def get_status_from_workflow(workflow_name, config_path, meta_hash):
    try:
        module_name = f"mindscape.dagtoy.workflows.{snake_case(workflow_name.removesuffix('Workflow'))}"
        module = importlib.import_module(module_name)
        klass = getattr(module, workflow_name)
        instance = klass(config_path=config_path, meta_hash=meta_hash)
        return instance.get_status()
    except Exception:
        return "not_started"

def render_dag_image(dag, meta_hash, config_path):
    dot = Digraph(comment="DAG Structure")
    project_path = Path(config_path).resolve().parents[0]
    logs_dir = project_path / "logs"
    color_map = {
        "completed": "green",
        "failed": "red",
        "in_progress": "yellow",
        "not_started": "gray"
    }

    # Detect subworkflow/root placeholder nodes (ending with "_subworkflow_placeholder" or "_subworkflow_root")
    subworkflow_parents = {}
    for node, deps in dag.items():
        for dep in deps:
            if dep.endswith("_subworkflow_placeholder") or dep.endswith("_subworkflow_root"):
                subworkflow_parents.setdefault(dep, set()).add(node)

    for node in dag:
        if node.endswith("_subworkflow_placeholder") or node.endswith("_subworkflow_root"):
            # Render subworkflow placeholders and roots in gray
            dot.node(node, style="filled", fillcolor="lightgray", fontcolor="black")
        else:
            status = get_status_from_workflow(node, config_path, meta_hash)
            dot.node(node, label=node, style="filled", fillcolor=color_map.get(status, "gray"))

    for node, deps in dag.items():
        for dep in deps:
            # If dependency is a subworkflow/root placeholder, draw edge from parent to placeholder/root
            if dep.endswith("_subworkflow_placeholder") or dep.endswith("_subworkflow_root"):
                dot.edge(node, dep, style="dashed", color="gray")
            else:
                dot.edge(dep, node)

    output_path = Path(config_path).parent / f"dag_{meta_hash}.png"
    dot.render(str(output_path.with_suffix("")), format="png", cleanup=True)
    print(f"🖼️ DAG visual with status saved to: {output_path}")

def run_pipeline(config_path):
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    # Parse dag from the workflows list in the config file
    dag = {}
    workflows = config.get("workflows", [])
    for wf in workflows:
        name = wf.get("name")
        deps = wf.get("depends_on", [])
        subworkflows = wf.get("subworkflows", [])
        if name:
            dag[name] = deps.copy()
            # Add a root placeholder node for subworkflow visibility, if any subworkflows are present
            if subworkflows:
                root_placeholder = f"{name}_subworkflow_root"
                dag[root_placeholder] = []
                dag[name].append(root_placeholder)
            # For each subworkflow, add a placeholder node and add as a child of the parent step
            for subwf in subworkflows:
                sub_node = f"{subwf}_subworkflow_placeholder"
                dag.setdefault(sub_node, [])
                dag[name].append(sub_node)
    meta_hash = compute_meta_hash(config_path)
    print(f"🔁 Meta hash for DAG run: {meta_hash}\n")

    render_dag_image(dag, meta_hash, config_path)

    order = topological_sort(dag)
    for name in order:
        if not name.endswith("Workflow"):
            raise ValueError(f"Workflow name '{name}' must end in 'Workflow'")

        base = name.removesuffix("Workflow")
        snake = snake_case(base)
        module = importlib.import_module(f"mindscape.dagtoy.workflows.{snake}")
        klass = getattr(module, name) 
        instance = klass(config_path=config_path, meta_hash=meta_hash)

        if not instance.is_complete():
            instance.run()
        else:
            print(f"✅ Skipping {name} (already complete)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_path", required=True, help="Path to the config.yaml")
    args = parser.parse_args()
    run_pipeline(args.config_path)