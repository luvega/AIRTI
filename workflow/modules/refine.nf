process REFINE {
    label 'gpu_boltz'

    input:
    path screened

    output:
    path 'boltz_top30.jsonl', emit: selected

    script:
    if (params.mock_tools) {
        """
        printf '%s\n' '{"ligand_id":"query-1","target_id":"P00533","boltz_score":0.86,"rank":1}' > boltz_top30.jsonl
        """
    } else {
        """
        airti-tf refine-boltz --candidates ${screened} --output boltz_top30.jsonl
        """
    }
}

