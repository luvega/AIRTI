process REPORT {
    label 'cpu_small'
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path simulated

    output:
    path 'final_report', emit: report_dir
    path 'job_status.sqlite', emit: state_db

    script:
    if (params.mock_tools) {
        """
        mkdir -p final_report
        python - ${simulated} <<'PY'
        import json
        import pathlib
        import sqlite3
        import sys

        records = [
            json.loads(line)
            for line in pathlib.Path(sys.argv[1]).read_text().splitlines()
            if line.strip()
        ]
        report_lines = [
            '# AI 反向钓靶编排模拟报告',
            '',
            '> 本报告仅验证 Nextflow 编排、批次身份传递和恢复，不包含真实计算结果。',
            '',
        ]
        for record in records:
            artifact_id = f"mock-md-{record['ligand_id']}"
            report_lines.append(
                f"- {record['ligand_id']}：候选 {record['target_id']}；"
                f"模拟 MD {record['completed_ns']} ns。 "
                f"<!--METRIC artifact:{artifact_id}-->"
            )
        report_path = pathlib.Path('final_report/report.md')
        report_path.write_text('\\n'.join(report_lines) + '\\n', encoding='utf-8')
        report_sha256 = __import__('hashlib').sha256(report_path.read_bytes()).hexdigest()
        manifest = {
            'schema_version': '1.0',
            'validation_scope': 'orchestration_mock_only',
            'query_count': len(records),
            'queries': [record['ligand_id'] for record in records],
            'technical_success_rate': 1.0 if records else 0.0,
            'all_metrics_traceable': True,
            'unsupported_targets_have_no_numeric_score': True,
            'report': {'path': 'report.md', 'sha256': report_sha256},
        }
        pathlib.Path('final_report/report_manifest.json').write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + '\\n',
            encoding='utf-8',
        )
        connection = sqlite3.connect('job_status.sqlite')
        connection.execute('CREATE TABLE tasks(ligand_id TEXT, stage TEXT, status TEXT)')
        connection.executemany(
            "INSERT INTO tasks VALUES (?, 'report', 'succeeded')",
            [(record['ligand_id'],) for record in records],
        )
        connection.commit()
        connection.close()
        PY
        """
    } else {
        """
        airti-tf render-report --candidates ${simulated} --output final_report --state-db job_status.sqlite
        """
    }
}
