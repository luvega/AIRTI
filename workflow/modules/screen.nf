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
        printf '%s\n' '{"ligand_id":"query-1","target_id":"P00533","calibrated_score":0.91,"rank":1}' > calibrated_top300.jsonl
        """
    } else {
        """
        airti-tf screen --ligands ${ligands} --targets ${targets} --output calibrated_top300.jsonl
        """
    }
}

