process SCREEN {
    label 'cpu_screen'

    input:
    path ligands
    path targets

    output:
    path 'calibrated_top300.jsonl', emit: calibrated

    script:
    if (params.mock_tools) {
        """
        python - ${ligands} <<'PY'
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
        pathlib.Path('calibrated_top300.jsonl').write_text(
            ''.join(json.dumps(row, sort_keys=True) + '\\n' for row in rows),
            encoding='utf-8',
        )
        PY
        """
    } else {
        """
        airti-tf screen --ligands ${ligands} --targets ${targets} --output calibrated_top300.jsonl
        """
    }
}
