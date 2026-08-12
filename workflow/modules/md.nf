process MD {
    label 'gpu_md'

    input:
    path refined

    output:
    path 'md_top10.jsonl', emit: completed

    script:
    if (params.mock_tools) {
        """
        printf '%s\n' '{"ligand_id":"query-1","target_id":"P00533","md_score":0.79,"completed_ns":100.0,"rank":1}' > md_top10.jsonl
        """
    } else {
        """
        airti-tf run-md --candidates ${refined} --output md_top10.jsonl --resume
        """
    }
}

