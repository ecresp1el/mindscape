from pathlib import Path

class DoubletScoringWorkflow:
    def __init__(self, project_path):
        self.project_path = Path(project_path)
        self.step_name = self.__class__.__name__
        self.output_dir = self.project_path / "results" / self.step_name
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self):
        print(f"🧪 Running {self.step_name}")
        (self.output_dir / "doublet_scoring_done.txt").write_text("Doublet scoring done.\n")