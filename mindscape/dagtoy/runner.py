import importlib
from dag_config import dag

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

def run_pipeline():
    order = topological_sort(dag)
    for name in order:
        if not name.endswith("Workflow"):
            raise ValueError(f"Workflow name '{name}' must end in 'Workflow'")

        base = name.removesuffix("Workflow")
        snake = snake_case(base)
        module = importlib.import_module(f"mindscape.dagtoy.workflows.{snake}")
        klass = getattr(module, name)
        instance = klass()

        if not instance.is_complete():
            instance.run()
        else:
            print(f"✅ Skipping {name} (already complete)")

if __name__ == "__main__":
    run_pipeline()
