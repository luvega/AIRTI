process MD {
    label 'gpu_md'

    input:
    path refined

    output:
    path 'md_top10.jsonl', emit: completed

    script:
    if (params.mock_tools) {
        """
        python - ${refined} <<'PY'
        import json
        import pathlib
        import sys

        rows = []
        for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
            if line.strip():
                row = json.loads(line)
                row.update({"md_score": 0.79, "completed_ns": 100.0, "mock_only": True})
                rows.append(row)
        pathlib.Path('md_top10.jsonl').write_text(
            ''.join(json.dumps(row, sort_keys=True) + '\\n' for row in rows),
            encoding='utf-8',
        )
        PY
        """
    } else {
        """
        airti-tf run-md --candidates ${refined} --output md_top10.jsonl --resume
        """
    }
}
