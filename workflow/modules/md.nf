process MD {
    label 'gpu_md'

    input:
    tuple path(refined), path(boltz_assets)

    output:
    tuple path('md_candidates.jsonl'), path('md_assets'), emit: completed

    script:
    if (params.mock_tools) {
        """
        mkdir -p md_assets
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
        pathlib.Path('md_candidates.jsonl').write_text(
            ''.join(json.dumps(row, sort_keys=True) + '\\n' for row in rows),
            encoding='utf-8',
        )
        PY
        """
    } else {
        """
        airti-tf run-md \
          --candidates ${refined} \
          --output md_candidates.jsonl \
          --asset-dir md_assets
        """
    }
}
