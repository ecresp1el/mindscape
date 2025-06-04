# mindscape/dagtoy/runner.py
import importlib
from dag_config import dag
import argparse
import hashlib

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

def run_pipeline(config_path):
    meta_hash = compute_meta_hash(config_path)
    print(f"🔁 Meta hash for DAG run: {meta_hash}\n")

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