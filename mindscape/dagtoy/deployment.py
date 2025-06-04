name: dev
flow_name: run_mindscape_flow
entrypoint: mindscape/dagtoy/prefect_flow.py:run_mindscape_flow
work_pool:
  name: default
  work_queue_name: default
tags: ["dagtoy"]