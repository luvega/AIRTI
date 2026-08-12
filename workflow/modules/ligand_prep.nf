process LIGAND_PREP {
    label 'cpu_small'

    input:
    path queries

    output:
    path 'prepared.smi', emit: prepared

    script:
    if (params.mock_tools) {
        """
        cp ${queries} prepared.smi
        count=\$(grep -cv '^[[:space:]]*\$' prepared.smi)
        test "\$count" -ge 1
        test "\$count" -le 5
        """
    } else {
        """
        airti-tf prepare-ligands ${queries} --output prepared.smi --max-molecules 5
        """
    }
}
