process TARGET_LIBRARY {
    label 'cpu_small'

    output:
    tuple path('targets.jsonl'), path('targets'), emit: library

    script:
    if (params.mock_tools) {
        """
        mkdir -p targets
        printf '%s\n' \
          '{"target_id":"P00533","status":"ready","input_sha256":"mock-ready"}' \
          '{"target_id":"P0UNSP","status":"unsupported","unsupported_reason":"no_structure","input_sha256":"mock-unsupported"}' \
          > targets.jsonl
        """
    } else {
        error "TARGET_LIBRARY is only used by the test profile"
    }
}
