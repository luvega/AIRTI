process LIGAND_PREP {
    label 'cpu_small'

    input:
    path queries

    output:
    tuple path('prepared_ligands.jsonl'), path('prepared_ligands'), emit: prepared

    script:
    if (params.mock_tools) {
        """
        mkdir -p prepared_ligands
        cp ${queries} prepared_ligands.jsonl
        count=\$(grep -cv '^[[:space:]]*\$' prepared_ligands.jsonl)
        test "\$count" -ge 1
        test "\$count" -le 5
        """
    } else {
        """
        airti-tf prepare-ligands ${queries} \
          --output prepared_ligands.jsonl \
          --asset-dir prepared_ligands \
          --max-molecules 5
        """
    }
}
