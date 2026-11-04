import os
import tempfile
import unittest

from final_report_utils import create_batch_index, maybe_create_batch_artifacts
from main import (
    _BATCH_SUMMARY_FIELDS,
    _run_command,
    _run_stage,
    _set_tool_status,
    _write_batch_issues_tsv,
    extract_batch_summary_fields,
)


class RunCommandErrorTests(unittest.TestCase):
    def test_failure_message_includes_stderr_tail(self):
        with self.assertRaises(RuntimeError) as ctx:
            _run_command(['sh', '-c', 'echo first >&2; echo boom >&2; exit 3'], 'FakeTool', sample_tag='s1')
        message = str(ctx.exception)
        self.assertIn('FakeTool failed for s1 with exit code 3', message)
        self.assertIn('boom', message)

    def test_missing_executable_is_reported(self):
        with self.assertRaises(RuntimeError) as ctx:
            _run_command(['afluid-no-such-tool-xyz'], 'Ghost')
        self.assertIn('executable not found', str(ctx.exception))

    def test_redirected_stderr_points_to_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            err_fp = os.path.join(tmpdir, 'tool.err')
            with open(err_fp, 'w') as err:
                with self.assertRaises(RuntimeError) as ctx:
                    _run_command(['sh', '-c', 'echo bad >&2; exit 1'], 'Logged', stderr=err)
            self.assertIn(err_fp, str(ctx.exception))


class StageStatusTests(unittest.TestCase):
    def test_stage_records_completed_and_failed(self):
        flags = {}
        self.assertEqual(_run_stage(flags, 'ok_stage', lambda: 42), 42)
        self.assertEqual(flags['Tool Status']['ok_stage']['status'], 'completed')

        persisted = []

        def boom():
            raise ValueError('kaput')

        with self.assertRaises(ValueError) as ctx:
            _run_stage(flags, 'c_blast', boom, on_failure=lambda: persisted.append(True))
        self.assertEqual(flags['Tool Status']['c_blast'], {'status': 'failed', 'detail': 'kaput'})
        self.assertEqual(ctx.exception.afluid_stage, 'c_blast')
        self.assertIs(ctx.exception.afluid_flags, flags)
        self.assertEqual(persisted, [True])


class BatchSummaryToolErrorTests(unittest.TestCase):
    def test_summary_distinguishes_failed_skipped_completed(self):
        flags = {'Sample': {}}
        _set_tool_status(flags, 'flumut', 'completed')
        _set_tool_status(flags, 'nextclade', 'failed', 'Nextclade H5 failed with exit code 1: dataset\tnot found')
        _set_tool_status(flags, 'genin', 'skipped', 'not routed')
        _set_tool_status(flags, 'flumut_db_update', 'failed', 'network down; retry later')

        summary = extract_batch_summary_fields(flags)

        self.assertEqual(summary['flumut_status'], 'completed')
        self.assertEqual(summary['nextclade_status'], 'failed')
        self.assertEqual(summary['genin2_status'], 'skipped')
        self.assertEqual(summary['getref_status'], '')
        self.assertEqual(summary['nextclade_ran'], 'no')
        self.assertIn('nextclade: Nextclade H5 failed with exit code 1: dataset not found', summary['tool_errors'])
        self.assertIn('flumut_db_update: network down, retry later', summary['tool_errors'])
        self.assertNotIn('\t', summary['tool_errors'])
        for field in ('flumut_status', 'nextclade_status', 'genin2_status', 'getref_status', 'tool_errors'):
            self.assertIn(field, _BATCH_SUMMARY_FIELDS)


class BatchIndexToolErrorTests(unittest.TestCase):
    def test_index_shows_tool_errors_and_batch_issues(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            row = {
                'input_file': 'sample1.fasta',
                'file_tag': 'sample1',
                'status': 'ok',
                'reports_dir_rel': 'sample1',
                'error': '',
                'tool_errors': 'nextclade: exit code 1; genin: boom',
            }
            index_fp = create_batch_index(
                tmpdir, [row], [row], 'index.html',
                batch_issues=[('batch_flumut_workbook', 'FluMut failed')],
            )
            with open(index_fp, encoding='utf-8') as handle:
                html = handle.read()

        self.assertIn('tool errors', html)
        self.assertIn('<li>nextclade: exit code 1</li>', html)
        self.assertIn('<li>genin: boom</li>', html)
        self.assertIn('Batch-level issues', html)
        self.assertIn('batch_flumut_workbook', html)
        self.assertIn('Successful with tool errors: 1', html)

    def test_index_is_written_when_every_sample_failed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_fp = os.path.join(tmpdir, 'batch_summary_allfail.tsv')
            with open(summary_fp, 'w', encoding='utf-8') as handle:
                handle.write('input_file\tfile_tag\tstatus\terror\treports_dir_rel\n')
                handle.write('s1.fasta\ts1\tfailed\tpreprocess: no sequences\ts1\n')

            created = maybe_create_batch_artifacts(tmpdir, summary_fp)

            self.assertTrue(created)
            index_fp = os.path.join(tmpdir, 'batch_index_allfail.html')
            self.assertTrue(os.path.isfile(index_fp))
            with open(index_fp, encoding='utf-8') as handle:
                self.assertIn('preprocess: no sequences', handle.read())
            self.assertFalse(os.path.exists(os.path.join(tmpdir, 'batch_reports_allfail.zip')))

    def test_batch_issues_tsv_flattens_whitespace(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fp = os.path.join(tmpdir, 'issues.tsv')
            _write_batch_issues_tsv(fp, [('batch_artifacts', 'line one\n\tline two')])
            with open(fp, encoding='utf-8') as handle:
                lines = handle.read().splitlines()
        self.assertEqual(lines, ['step\terror', 'batch_artifacts\tline one line two'])


if __name__ == '__main__':
    unittest.main()
