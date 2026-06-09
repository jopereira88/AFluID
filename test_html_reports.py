import os
import tempfile
import unittest

from final_report_utils import (
    apply_tool_version_placeholder,
    create_batch_index,
    create_flumut_table_multi,
    create_sample_table,
    create_segment_background_table,
    create_segment_background_table_multi,
    create_segment_mutations_table,
    html_skeleton,
    html_skeleton_multi,
)
from metadata_tree_utils import load_reference_dicts
from structures import seg_lens


class HtmlReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geo_dict, cls.taxa_dict = load_reference_dicts(
            'metadata',
            'geo_tree.json',
            'host_taxonomy_tree.json',
        )

    def test_batch_index_moves_id_report_and_flags_to_extras(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, 'sample1')
            os.mkdir(sample_dir)

            for name in (
                'sample1_final_report.html',
                'sample1_ID_Report.txt',
                'sample1_flags.json',
                'sample1_extra.tsv',
            ):
                with open(os.path.join(sample_dir, name), 'w', encoding='utf-8') as handle:
                    handle.write('x')

            row = {
                'input_file': 'sample1.fasta',
                'file_tag': 'sample1',
                'status': 'ok',
                'reports_dir_rel': 'sample1',
                'error': '',
            }
            index_fp = create_batch_index(tmpdir, [row], [row], 'index.html')

            with open(index_fp, 'r', encoding='utf-8') as handle:
                html = handle.read()

            self.assertNotIn('<th>ID report</th>', html)
            self.assertNotIn('<th>Flags</th>', html)
            self.assertIn('sample1_ID_Report.txt', html)
            self.assertIn('sample1_flags.json', html)
            self.assertIn('sample1_extra.tsv', html)
            self.assertIn('target="_blank"', html)
            self.assertIn('rel="noopener noreferrer"', html)

    def test_sample_table_uses_canonical_segment_order_and_locus_mutations(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            id_report_fp = os.path.join(tmpdir, 'sample_ID_Report.txt')
            with open(id_report_fp, 'w', encoding='utf-8') as handle:
                handle.write('SAMPLE_NAME\n')
                handle.write('>ha_name\n')
                handle.write('>pb1_name\n')

            flags = {
                'Sample': {
                    'PB2': [],
                    'PB1': ['>int_pb1'],
                    'PA': [],
                    'HA': ['>int_ha'],
                    'NP': [],
                    'NA': [],
                    'MP': [],
                    'NS': [],
                    'PB1_len': 2220,
                    'HA_len': 1700,
                    'PB1_ref': ['CY000002.1'],
                    'HA_ref': ['CY000004.1'],
                    'PB1_muts': ['13P'],
                    'HA_muts': ['156A'],
                }
            }
            mappings = {'>int_pb1': '>pb1_name', '>int_ha': '>ha_name'}

            html = create_sample_table(id_report_fp, flags, seg_lens, mappings, {'13P': ('PB1:13P', ''), '156A': ('HA:156A', '')})

            ordered_segments = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
            positions = [html.find(f'<td>{segment}</td>') for segment in ordered_segments]
            self.assertTrue(all(position >= 0 for position in positions))
            self.assertEqual(positions, sorted(positions))
            self.assertIn('PB1:13P', html)
            self.assertIn('HA:156A', html)
            self.assertIn('https://www.ncbi.nlm.nih.gov/nuccore/CY000002.1', html)
            self.assertIn('https://www.ncbi.nlm.nih.gov/nuccore/CY000004.1', html)

    def test_single_segment_background_sorts_observed_rows_and_links_references(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            id_report_fp = os.path.join(tmpdir, 'sample_ID_Report.txt')
            with open(id_report_fp, 'w', encoding='utf-8') as handle:
                handle.write('SAMPLE_NAME\tSEGMENT\tCLUSTER\tREPRESENTATIVE\t%ID\tASSIGNED_BY\tGENOTYPE\tHOST\tCOUNTRY\tCOLLECTION_DATE\n')
                handle.write('>ha_name\t4\tUNCLUSTERED\tCY000004.1\t98\tL-BLAST\tH1N1\tHomo sapiens\tBrazil\t2024-02\n')
                handle.write('>pb1_name\t2\tUNCLUSTERED\tCY000002.1\t99\tL-BLAST\tH1N1\tHomo sapiens\tBrazil\t2024-01\n')

            flags = {
                'Sample': {
                    'PB2': [],
                    'PB1': ['>int_pb1'],
                    'PA': [],
                    'HA': ['>int_ha'],
                    'NP': [],
                    'NA': [],
                    'MP': [],
                    'NS': [],
                    'PB1_len': 2220,
                    'HA_len': 1700,
                    'PB1_ref': ['CY000002.1'],
                    'HA_ref': ['CY000004.1'],
                }
            }
            mappings = {'>int_pb1': '>pb1_name', '>int_ha': '>ha_name'}

            html = create_segment_background_table(
                id_report_fp,
                flags,
                seg_lens,
                mappings,
                self.geo_dict,
                self.taxa_dict,
            )

            self.assertLess(html.find('<td>PB1</td>'), html.find('<td>HA</td>'))
            self.assertLess(html.find('<th>Note</th>'), html.find('<th>Cluster</th>'))
            self.assertLess(html.find('<th>Cluster</th>'), html.find('<th>Reference</th>'))
            self.assertLess(html.find('<th>Assigned_By</th>'), html.find('<th>Genotype</th>'))
            self.assertIn('https://www.ncbi.nlm.nih.gov/nuccore/CY000002.1', html)
            self.assertIn('https://www.ncbi.nlm.nih.gov/nuccore/CY000004.1', html)

    def test_multi_segment_background_sorts_by_canonical_segment_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            id_report_fp = os.path.join(tmpdir, 'sample_ID_Report.txt')
            formatted_fasta_fp = os.path.join(tmpdir, 'format_sample.fasta')

            with open(id_report_fp, 'w', encoding='utf-8') as handle:
                handle.write('SAMPLE_NAME\tSEGMENT\tCLUSTER\tREPRESENTATIVE\t%ID\tASSIGNED_BY\tGENOTYPE\tHOST\tCOUNTRY\tCOLLECTION_DATE\n')
                handle.write('ha_seq\t4\tUNCLUSTERED\tCY000004.1\t98\tL-BLAST\tH1N1\tHomo sapiens\tBrazil\t2024-02\n')
                handle.write('pb1_seq\t2\tUNCLUSTERED\tCY000002.1\t99\tL-BLAST\tH1N1\tHomo sapiens\tBrazil\t2024-01\n')

            with open(formatted_fasta_fp, 'w', encoding='utf-8') as handle:
                handle.write('>int_ha\nAAAA\n>int_pb1\nAAAA\n')

            mappings = {'>int_ha': 'ha_seq', '>int_pb1': 'pb1_seq'}

            html = create_segment_background_table_multi(
                id_report_fp,
                seg_lens,
                mappings,
                formatted_fasta_fp,
                self.geo_dict,
                self.taxa_dict,
            )

            self.assertLess(html.find('<td>PB1</td>'), html.find('<td>HA</td>'))
            self.assertLess(html.find('<th>SeqLen</th>'), html.find('<th>Note</th>'))
            self.assertLess(html.find('<th>Note</th>'), html.find('<th>Cluster</th>'))
            self.assertLess(html.find('<th>Assigned_By</th>'), html.find('<th>Genotype</th>'))
            self.assertIn('https://www.ncbi.nlm.nih.gov/nuccore/CY000002.1', html)
            self.assertIn('target="_blank"', html)

    def test_apply_tool_version_placeholder_replaces_custom_tag(self):
        html = '<div><!--VERSIONS--></div>'
        updated = apply_tool_version_placeholder(
            html,
            tool_versions={'cd-hit': {'version': '4.8.1', 'build': 'Apr 24 2025'}},
            placeholder_tag='<!--VERSIONS-->',
        )

        self.assertIn('cd-hit: 4.8.1 Apr 24 2025', updated)
        self.assertNotIn('<!--VERSIONS-->', updated)

    def test_apply_tool_version_placeholder_replaces_per_tool_tags(self):
        html = (
            '<div>'
            '<span><!--TOOL_VERSION_BLASTN--></span>'
            '<span><!--TOOL_VERSION_CD_HIT_EST_2D--></span>'
            '<span><!--TOOL_VERSION_FLUMUT--></span>'
            '<span><!--TOOL_VERSION_FLUMUTDB--></span>'
            '<span><!--TOOL_VERSION_GENIN2--></span>'
            '<span><!--TOOL_VERSION_NEXTCLADE--></span>'
            '</div>'
        )
        updated = apply_tool_version_placeholder(
            html,
            tool_versions={
                'blastn': '2.16.0+',
                'cd-hit-est-2d': '4.8.1 (built on Apr 24 2025)',
                'flumut': '1.2.3',
                'FluMutDB': '2026.05',
                'genin2': '0.7.0',
                'nextclade': '3.12.0',
            },
        )

        self.assertIn('<span>2.16.0+</span>', updated)
        self.assertIn('<span>4.8.1 (built on Apr 24 2025)</span>', updated)
        self.assertIn('<span>1.2.3</span>', updated)
        self.assertIn('<span>2026.05</span>', updated)
        self.assertIn('<span>0.7.0</span>', updated)
        self.assertIn('<span>3.12.0</span>', updated)
        self.assertNotIn('<!--TOOL_VERSION_BLASTN-->', updated)

    def test_apply_tool_version_placeholder_fills_missing_per_tool_tags_with_unknown(self):
        html = '<div><!--TOOL_VERSION_BLASTN-->|<!--TOOL_VERSION_NEXTCLADE--></div>'
        updated = apply_tool_version_placeholder(html, tool_versions={'blastn': '2.16.0+'})

        self.assertIn('<div>2.16.0+|Unknown</div>', updated)

    def test_segment_mutations_table_uses_requested_locus_order(self):
        flags = {
            'Sample': {
                'PB2_muts': ['627K'],
                'PB1_muts': ['13P', '66S'],
                'PA_muts': ['38M', '127V'],
                'HA_muts': ['156A'],
                'NP_muts': ['105V'],
                'NA_muts': ['275Y'],
                'MP_muts': ['31N', '95K'],
                'NS_muts': [],
            }
        }
        lookup = {
            '627K': ('PB2:627K', ''),
            '13P': ('PB1:13P', ''),
            '66S': ('PB1-F2:66S', ''),
            '38M': ('PA:38M', ''),
            '127V': ('PA-X:127V', ''),
            '156A': ('HA:156A', ''),
            '105V': ('NP:105V', ''),
            '275Y': ('NA:275Y', ''),
            '95K': ('M1:95K', ''),
            '31N': ('M2:31N', ''),
        }

        html = create_segment_mutations_table(flags, lookup)

        ordered_labels = [
            'PB2:627K',
            'PB1:13P',
            'PB1-F2:66S',
            'PA:38M',
            'PA-X:127V',
            'HA:156A',
            'NP:105V',
            'NA:275Y',
            'M1:95K',
            'M2:31N',
        ]
        positions = [html.find(f'<td>{label}</td>') for label in ordered_labels]
        self.assertTrue(all(position >= 0 for position in positions))
        self.assertEqual(positions, sorted(positions))

    def test_flumut_multi_table_groups_by_sample_then_requested_locus_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            markers_fp = os.path.join(tmpdir, 'batch_markers.tsv')
            with open(markers_fp, 'w', encoding='utf-8') as handle:
                handle.write('Sample\tMutations in your sample\n')
                handle.write('sample_b\tM2:A31N\n')
                handle.write('sample_a\tM1:T95K\n')
                handle.write('sample_a\tPB1-F2:N66S\n')
                handle.write('sample_a\tPB2:E627K\n')
                handle.write('sample_b\tPB1:I13P\n')

            html = create_flumut_table_multi(
                reports_p=tmpdir,
                file_tag='batch',
                muts_interest={
                    'PB2': ['627K'],
                    'PB1': ['13P', '66S'],
                    'MP': ['95K', '31N'],
                },
                muts_loci_meaning={
                    '627K': ('PB2:627K', ''),
                    '13P': ('PB1:13P', ''),
                    '66S': ('PB1-F2:66S', ''),
                    '95K': ('M1:95K', ''),
                    '31N': ('M2:31N', ''),
                },
            )

            sample_a_pb2 = html.find('<td>sample_a</td>')
            sample_a_pb1f2 = html.find('<td>PB1-F2</td>')
            sample_a_m1 = html.find('<td>M1</td>')
            sample_b = html.find('<td>sample_b</td>')
            sample_b_pb1 = html.find('<td>PB1</td>', sample_b)
            sample_b_m2 = html.find('<td>M2</td>', sample_b)

            self.assertTrue(all(position >= 0 for position in (sample_a_pb2, sample_a_pb1f2, sample_a_m1, sample_b, sample_b_pb1, sample_b_m2)))
            self.assertLess(sample_a_pb2, sample_a_pb1f2)
            self.assertLess(sample_a_pb1f2, sample_a_m1)
            self.assertLess(sample_a_m1, sample_b)
            self.assertLess(sample_b_pb1, sample_b_m2)

    def test_single_sample_skeleton_uses_supported_tool_version_tags(self):
        self.assertIn('<!--TOOL_VERSION_BLASTN-->', html_skeleton)
        self.assertIn('<!--TOOL_VERSION_NEXTCLADE-->', html_skeleton)
        self.assertIn('<!--TOOL_VERSION_FLUMUT-->', html_skeleton)
        self.assertIn('<!--TOOL_VERSION_FLUMUTDB-->', html_skeleton)
        self.assertNotIn('<!--TOOL_VERSION_BLAST-->', html_skeleton)
        self.assertNotIn('<!--TOOL_VERSION_NETCLADE-->', html_skeleton)
        self.assertIn('<h2>Sample Information</h2>', html_skeleton)
        self.assertIn('<h2>Segment Background</h2>', html_skeleton)
        self.assertIn('<h2>Segment Mutations</h2>', html_skeleton)

    def test_multi_sample_skeleton_uses_supported_tool_version_tags(self):
        self.assertIn('<!--TOOL_VERSION_CD_HIT_EST_2D-->', html_skeleton_multi)
        self.assertIn('<!--TOOL_VERSION_BLASTN-->', html_skeleton_multi)
        self.assertIn('<!--TOOL_VERSION_NEXTCLADE-->', html_skeleton_multi)
        self.assertIn('<!--TOOL_VERSION_FLUMUT-->', html_skeleton_multi)
        self.assertIn('<!--TOOL_VERSION_FLUMUTDB-->', html_skeleton_multi)

    def test_flumut_sections_link_to_github_repository_in_new_tab(self):
        for html in (html_skeleton, html_skeleton_multi):
            self.assertIn('https://github.com/jopereira88/AFluID', html)
            self.assertIn('AFluID GitHub repository', html)
            self.assertIn('target="_blank"', html)
            self.assertIn('rel="noopener noreferrer"', html)

    def test_genin2_panel_text_uses_version_tag_in_single_and_multi_renderers(self):
        single_html = html_skeleton
        single_html = single_html.replace('<!--TAB4_INPUT-->', '<input type="radio" id="tab4" name="tab">')
        single_html = single_html.replace('<!--TAB4_LABEL-->', '<label for="tab4">Genin2</label>')
        single_html = single_html.replace('<!--TAB4_ACTIVE_CSS-->', '')
        single_html = single_html.replace('<!--TAB4_DISPLAY_CSS-->', '')
        single_html = single_html.replace(
            '<!--TAB4_PANEL-->',
            '''
        <article id="panel4" class="tab-panel">
          <h2>Genin2</h2>
          <p class="muted">
            GenIn2 performs genotype constellation assignment. This tab reports the
            constellation assignment for this sample using GenIn2 (<!--TOOL_VERSION_GENIN2-->).
          </p>
          <!--GENIN2_TABLE-->
        </article>
        '''
        )
        self.assertIn('<!--TOOL_VERSION_GENIN2-->', single_html)

        multi_html = html_skeleton_multi.replace(
            '<!--TAB4_PANEL-->',
            '''
        <article id="panel4" class="tab-panel">
          <h2>Genin2</h2>
          <p class="muted">
            GenIn2 performs genotype constellation assignment. This tab reports
            constellation assignments for eligible samples in this analysis using
            GenIn2 (<!--TOOL_VERSION_GENIN2-->).
          </p>
          <div>table</div>
        </article>
        '''
        )
        self.assertIn('<!--TOOL_VERSION_GENIN2-->', multi_html)


if __name__ == '__main__':
    unittest.main()
