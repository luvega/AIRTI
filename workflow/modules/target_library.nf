process TARGET_LIBRARY {
    label 'cpu_small'

    output:
    path 'targets.jsonl', emit: manifest

    script:
    if (params.mock_tools) {
        """
        printf '%s\n' \
          '{"target_id":"P00533","status":"ready","input_sha256":"mock-ready"}' \
          '{"target_id":"P0UNSP","status":"unsupported","unsupported_reason":"no_structure","input_sha256":"mock-unsupported"}' \
          > targets.jsonl
        """
    } else {
        """
        test -s ${params.target_manifest}
        cp ${params.target_manifest} targets.jsonl
        """
    }
}
