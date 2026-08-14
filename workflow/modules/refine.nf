process REFINE {
    label 'gpu_boltz'

    input:
    tuple path(screened), path(screen_assets)

    output:
    tuple path('boltz_candidates.jsonl'), path('boltz_assets'), emit: selected

    script:
    if (params.mock_tools) {
        """
        mkdir -p boltz_assets
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
        pathlib.Path('boltz_candidates.jsonl').write_text(
            ''.join(json.dumps(row, sort_keys=True) + '\\n' for row in rows),
            encoding='utf-8',
        )
        PY
        """
    } else {
        """
        airti-tf refine-boltz \
          --candidates ${screened} \
          --output boltz_candidates.jsonl \
          --asset-dir boltz_assets \
          --profile production \
          --top-n ${params.boltz_top_n} \
          --cache /models/boltz
        """
    }
}
