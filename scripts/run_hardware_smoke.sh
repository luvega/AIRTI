#!/usr/bin/env bash
set -euo pipefail

IMAGE="${AIRTI_IMAGE:-airti-tf:0.1.0-gpu}"
CACHE_ROOT="${AIRTI_BOLTZ_CACHE:-/mnt/ssd4t/airti-target-fishing/boltz}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
FIXTURE_ROOT="${REPO_ROOT}/tests/fixtures/smoke"
RUN_ROOT="${1:-${REPO_ROOT}/results/hardware-smoke/$(date +%Y%m%d-%H%M%S)}"
HOST_USER="$(id -u):$(id -g)"

CONFIDENCE_SHA="090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1"
AFFINITY_SHA="dcc5cd3722b1c9eaa34267e4ae32f55cbbf1963f4c19319381ccfa30fdd2ca9e"

if [[ -e "${RUN_ROOT}" ]]; then
    echo "Refusing to overwrite existing smoke directory: ${RUN_ROOT}" >&2
    exit 64
fi

mkdir -p "${RUN_ROOT}/docking" "${RUN_ROOT}/gromacs" "${RUN_ROOT}/boltz"
mkdir -p "${CACHE_ROOT}"

docker image inspect "${IMAGE}" --format '{{.Id}}' > "${RUN_ROOT}/image-digest.txt"

docker run --rm \
    --user "${HOST_USER}" \
    -v "${FIXTURE_ROOT}/docking:/input:ro" \
    -v "${RUN_ROOT}/docking:/output" \
    "${IMAGE}" \
    qvina2 --receptor /input/receptor.pdbqt --ligand /input/ligand.pdbqt \
    --center_x 0 --center_y 0 --center_z 0 \
    --size_x 12 --size_y 12 --size_z 12 \
    --exhaustiveness 1 --num_modes 1 --seed 29 \
    --out /output/qvina-out.pdbqt --log /output/qvina.log \
    2>&1 | tee "${RUN_ROOT}/docking/qvina.stdout.log"
test -s "${RUN_ROOT}/docking/qvina-out.pdbqt"
grep -Fq "Writing output ... done." "${RUN_ROOT}/docking/qvina.stdout.log"

docker run --rm \
    --user "${HOST_USER}" \
    -v "${FIXTURE_ROOT}/gromacs:/input:ro" \
    -v "${RUN_ROOT}/gromacs:/output" \
    "${IMAGE}" sh -c \
    'cd /output && cp /input/topol.top topol.top && gmx solvate -cs spc216.gro -box 2 2 2 -o solvated.gro -p topol.top && gmx grompp -f /input/md.mdp -c solvated.gro -p topol.top -o md.tpr -maxwarn 1' \
    2>&1 | tee "${RUN_ROOT}/gromacs/setup.log"

docker run --rm --gpus all \
    --user "${HOST_USER}" \
    -v "${RUN_ROOT}/gromacs:/output" \
    "${IMAGE}" sh -c \
    'cd /output && gmx mdrun -s md.tpr -nsteps 1000 -deffnm md -nb gpu -pin on -ntmpi 1 -ntomp 4' \
    2>&1 | tee "${RUN_ROOT}/gromacs/mdrun.stdout.log"
grep -Fq "1 GPU selected for this run." "${RUN_ROOT}/gromacs/mdrun.stdout.log"
grep -Fq "Finished mdrun" "${RUN_ROOT}/gromacs/md.log"

docker run --rm --gpus all \
    --user "${HOST_USER}" \
    -e HOME=/tmp/airti-home \
    -v "${FIXTURE_ROOT}/boltz:/input:ro" \
    -v "${CACHE_ROOT}:/models/boltz" \
    -v "${RUN_ROOT}/boltz:/output" \
    "${IMAGE}" \
    boltz predict /input/input.yaml --out_dir /output --cache /models/boltz \
    --model boltz2 --diffusion_samples 1 --sampling_steps 10 \
    --diffusion_samples_affinity 1 --sampling_steps_affinity 10 \
    --max_parallel_samples 1 --num_workers 0 --seed 29 \
    2>&1 | tee "${RUN_ROOT}/boltz/boltz.stdout.log"
grep -Fq "GPU available: True" "${RUN_ROOT}/boltz/boltz.stdout.log"
grep -Fq "Number of failed examples: 0" "${RUN_ROOT}/boltz/boltz.stdout.log"
test -s "${RUN_ROOT}/boltz/boltz_results_input/predictions/input/input_model_0.cif"
test -s "${RUN_ROOT}/boltz/boltz_results_input/predictions/input/affinity_input.json"

printf '%s  %s\n' "${CONFIDENCE_SHA}" "${CACHE_ROOT}/boltz2_conf.ckpt" | sha256sum -c -
printf '%s  %s\n' "${AFFINITY_SHA}" "${CACHE_ROOT}/boltz2_aff.ckpt" | sha256sum -c -

find "${RUN_ROOT}" -type f ! -name artifacts.sha256 -print0 \
    | sort -z \
    | xargs -0 sha256sum > "${RUN_ROOT}/artifacts.sha256"

echo "AIRTI hardware smoke passed: ${RUN_ROOT}"
