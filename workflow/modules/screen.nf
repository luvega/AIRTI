process SCREEN {
    label 'cpu_screen'

    input:
    tuple path(ligand_manifest), path(ligand_assets)
    tuple path(target_manifest), path(target_assets)

    output:
    tuple path('screened_candidates.jsonl'), path('screen_assets'), emit: calibrated

    script:
    if (params.mock_tools) {
        """
        mkdir -p screen_assets
        python - ${ligand_manifest} <<'PY'
        import json
        import pathlib
        import sys

        rows = []
        for index, line in enumerate(pathlib.Path(sys.argv[1]).read_text().splitlines(), 1):
            if not line.strip():
                continue
            fields = line.split()
            ligand_id = fields[1] if len(fields) > 1 else f"query-{index}"
            rows.append({
                "ligand_id": ligand_id,
                "target_id": "P00533",
                "calibrated_score": 0.91,
                "rank": 1,
                "mock_only": True,
            })
        pathlib.Path('screened_candidates.jsonl').write_text(
            ''.join(json.dumps(row, sort_keys=True) + '\\n' for row in rows),
            encoding='utf-8',
        )
        PY
        """
    } else {
        """
        airti-tf screen \
          --ligands ${ligand_manifest} \
          --targets ${target_manifest} \
          --output screened_candidates.jsonl \
          --asset-dir screen_assets \
          --top-n ${params.screen_top_n}
        """
    }
}
