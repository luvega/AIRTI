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

        record = json.loads(pathlib.Path(sys.argv[1]).read_text().splitlines()[0])
        pathlib.Path('final_report/report.md').write_text(
            '# AI 反向钓靶模拟报告\\n\\n'
            f"候选靶点：{record['target_id']}；MD 完成：{record['completed_ns']} ns。\\n",
            encoding='utf-8',
        )
        connection = sqlite3.connect('job_status.sqlite')
        connection.execute('CREATE TABLE tasks(stage TEXT, status TEXT)')
        connection.execute("INSERT INTO tasks VALUES ('report', 'succeeded')")
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

