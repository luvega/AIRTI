process EVALUATE_CASE {
    label 'cpu_small'
    publishDir "${params.outdir}/case_evaluation", mode: 'copy', overwrite: true

    input:
    tuple path(screen_manifest), path(screen_assets)
    tuple path(boltz_manifest), path(boltz_assets)
    tuple path(md_manifest), path(md_assets)
    path case_definition

    output:
    path 'case_evaluation.json', emit: evaluation

    script:
    if (params.mock_tools) {
        """
        printf '%s\n' '{"schema_version":"1.0","validation_scope":"orchestration_mock_only"}' \
          > case_evaluation.json
        """
    } else {
        """
        airti-tf evaluate-case \
          --case ${case_definition} \
          --screen ${screen_manifest} \
          --boltz ${boltz_manifest} \
          --md ${md_manifest} \
          --output case_evaluation.json
        """
    }
}
