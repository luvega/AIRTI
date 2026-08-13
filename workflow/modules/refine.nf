process REFINE {
    label 'gpu_boltz'

    input:
    path screened

    output:
    path 'boltz_top30.jsonl', emit: selected

    script:
    if (params.mock_tools) {
        """
        python - ${screened} <<'PY'
        import json
        import pathlib
        import sys

        rows = []
        for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
            if line.strip():
                row = json.loads(line)
                row.update({"boltz_score": 0.86, "mock_only": True})
                rows.append(row)
        pathlib.Path('boltz_top30.jsonl').write_text(
            ''.join(json.dumps(row, sort_keys=True) + '\\n' for row in rows),
            encoding='utf-8',
        )
        PY
        """
    } else {
        """
        airti-tf refine-boltz --candidates ${screened} --output boltz_top30.jsonl
        """
    }
}
