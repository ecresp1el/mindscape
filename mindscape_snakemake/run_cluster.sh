#!/bin/bash
#SBATCH --job-name={{rule}}_{{wildcards}}
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/MindscapeProjects/logs/{{rule}}_{{wildcards}}_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/MindscapeProjects/logs/{{rule}}_{{wildcards}}_%j.err
#SBATCH --time={{resources.time}}
#SBATCH --mem={{resources.mem_mb}}
#SBATCH --cpus-per-task={{resources.cpus}}
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elcrespo@umich.edu

set -euo pipefail
echo "🚀 Running rule: {{rule}} | Wildcards: {{wildcards}} | Host: $HOSTNAME"
{{exec_job}}