nextflow.enable.dsl = 2

include { TARGET_LIBRARY } from './modules/target_library'
include { LIGAND_PREP } from './modules/ligand_prep'
include { SCREEN } from './modules/screen'
include { REFINE } from './modules/refine'
include { MD } from './modules/md'
include { REPORT } from './modules/report'

workflow {
    if (!params.queries) {
        error "--queries is required (SDF or SMILES; 1-5 molecules)"
    }

    query_ch = Channel.fromPath(params.queries, checkIfExists: true)
    if (params.mock_tools) {
        targets = TARGET_LIBRARY().library
    } else {
        target_manifest_ch = Channel.fromPath(
            params.target_manifest,
            checkIfExists: true,
        )
        target_assets_ch = Channel.fromPath(
            params.target_assets,
            checkIfExists: true,
        )
        targets = target_manifest_ch.combine(target_assets_ch)
    }
    ligands = LIGAND_PREP(query_ch)
    screened = SCREEN(ligands.prepared, targets)
    refined = REFINE(screened.calibrated)
    simulated = MD(refined.selected)
    REPORT(simulated.completed)
}
