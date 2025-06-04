from prefect import task
from mindscape.dagtoy.workflows.data_import import DataImportWorkflow

@task
def run_data_import(config_path, meta_hash=None):
    print("🔧 Running DataImportWorkflow via Prefect")
    wf = DataImportWorkflow(config_path=config_path, meta_hash=meta_hash)
    if not wf.is_complete():
        wf.run()