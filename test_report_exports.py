import os
import tempfile
import unittest

import pandas as pd

from final_report_utils import export_cluster_composition, export_id_report_rollup
from metadata_tree_utils import load_reference_dicts


class ReportExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geo_dict, cls.taxa_dict = load_reference_dicts(
            'metadata',
            'geo_tree.json',
            'host_taxonomy_tree.json',
        )

    def test_export_id_report_rollup_adds_static_rollup_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            id_report_fp = os.path.join(tmpdir, 'sample_ID_Report.txt')
            output_fp = os.path.join(tmpdir, 'sample_ID_Report_rollup.tsv')

            with open(id_report_fp, 'w', encoding='utf-8') as handle:
                handle.write(
                    'SAMPLE_NAME\tSEGMENT\tCLUSTER\tREPRESENTATIVE\t%ID\tASSIGNED_BY\tGENOTYPE\tHOST\tCOUNTRY\tCOLLECTION_DATE\n'
                )
                handle.write(
                    ">sample1\t{'4': 10}\t>Cluster 10\tREP1\t99.1\tCD-HIT\t{'H5N1': 1.0}\t{'Gallus gallus': 0.8, 'Anas platyrhynchos': 0.2}\t{'USA': 0.9, 'Brazil': 0.1}\t2024-01-02\n"
                )
                handle.write(
                    '>sample2\t4\tUNCLUSTERED\tREP2\t98.2\tL-BLAST\tH1N1\tHomo sapiens\tBrazil\t2024-02\n'
                )

            export_id_report_rollup(id_report_fp, output_fp, self.geo_dict, self.taxa_dict)

            df = pd.read_table(output_fp, index_col=False)

            self.assertEqual(
                list(df.columns[-4:]),
                ['HOST_FAMILY', 'HOST_CLASS', 'COUNTRY_REGION', 'COUNTRY_CONTINENT'],
            )
            self.assertEqual(df.loc[0, 'HOST_FAMILY'], "{'Phasianidae': 0.80, 'Anatidae': 0.20}")
            self.assertEqual(df.loc[0, 'HOST_CLASS'], "{'Aves': 1.00}")
            self.assertEqual(df.loc[0, 'COUNTRY_REGION'], "{'Northern America': 0.90, 'South America': 0.10}")
            self.assertEqual(df.loc[0, 'COUNTRY_CONTINENT'], "{'North America': 0.90, 'South America': 0.10}")
            self.assertEqual(df.loc[1, 'HOST_FAMILY'], 'Hominidae')
            self.assertEqual(df.loc[1, 'HOST_CLASS'], 'Mammalia')
            self.assertEqual(df.loc[1, 'COUNTRY_REGION'], 'South America')
            self.assertEqual(df.loc[1, 'COUNTRY_CONTINENT'], 'South America')

    def test_export_id_report_rollup_rounds_to_two_decimals_and_sums_to_one(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            id_report_fp = os.path.join(tmpdir, 'sample_ID_Report.txt')
            output_fp = os.path.join(tmpdir, 'sample_ID_Report_rollup.tsv')

            with open(id_report_fp, 'w', encoding='utf-8') as handle:
                handle.write(
                    'SAMPLE_NAME\tSEGMENT\tCLUSTER\tREPRESENTATIVE\t%ID\tASSIGNED_BY\tGENOTYPE\tHOST\tCOUNTRY\tCOLLECTION_DATE\n'
                )
                handle.write(
                    ">sample1\t{'4': 10}\t>Cluster 10\tREP1\t99.1\tCD-HIT\t{'H5N1': 1.0}\t{'Gallus gallus': 1, 'Anas platyrhynchos': 1, 'Homo sapiens': 1}\t{'USA': 1, 'France': 1, 'Brazil': 1}\t2024-01-02\n"
                )

            export_id_report_rollup(id_report_fp, output_fp, self.geo_dict, self.taxa_dict)

            df = pd.read_table(output_fp, index_col=False)

            self.assertEqual(df.loc[0, 'HOST_FAMILY'], "{'Anatidae': 0.34, 'Hominidae': 0.33, 'Phasianidae': 0.33}")
            self.assertEqual(df.loc[0, 'HOST_CLASS'], "{'Aves': 0.67, 'Mammalia': 0.33}")
            self.assertEqual(df.loc[0, 'COUNTRY_REGION'], "{'Northern America': 0.34, 'South America': 0.33, 'Western Europe': 0.33}")
            self.assertEqual(df.loc[0, 'COUNTRY_CONTINENT'], "{'Europe': 0.34, 'North America': 0.33, 'South America': 0.33}")

    def test_export_cluster_composition_uses_only_clusters_present_in_analysis(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_fp = os.path.join(tmpdir, 'metadata.csv')
            clstr_fp = os.path.join(tmpdir, 'clusters.clstr')
            id_report_fp = os.path.join(tmpdir, 'sample_ID_Report.txt')
            output_fp = os.path.join(tmpdir, 'sample_cluster_composition.tsv')

            with open(metadata_fp, 'w', encoding='utf-8') as handle:
                handle.write('ACCESSION;GENOTYPE;SEGMENT;HOST;COUNTRY;COLLECTION_DATE\n')
                handle.write('ACC1;H5N1;4;Gallus gallus;Germany;2024-01-01\n')
                handle.write('ACC2;H5N1;4;Anas platyrhynchos;France;2024-01-02\n')
                handle.write('ACC3;H1N1;6;Homo sapiens;Brazil;2024-01-03\n')

            with open(clstr_fp, 'w', encoding='utf-8') as handle:
                handle.write('>Cluster 1\n')
                handle.write('0\t100nt, >ACC1... *\n')
                handle.write('1\t100nt, >ACC2... at 99%\n')
                handle.write('>Cluster 2\n')
                handle.write('0\t100nt, >ACC3... *\n')

            with open(id_report_fp, 'w', encoding='utf-8') as handle:
                handle.write('SAMPLE_NAME\tSEGMENT\tCLUSTER\tREPRESENTATIVE\t%ID\tASSIGNED_BY\tGENOTYPE\tHOST\tCOUNTRY\n')
                handle.write(
                    ">sample1\t{'4': 2}\t>Cluster 1\tREP1\t99.1\tCD-HIT\t{'H5N1': 1.0}\t{'Gallus gallus': 0.5, 'Anas platyrhynchos': 0.5}\t{'Germany': 0.5, 'France': 0.5}\n"
                )
                handle.write('>sample2\t4\tUNCLUSTERED\tREP2\t98.0\tL-BLAST\tH1N1\tHomo sapiens\tBrazil\n')

            export_cluster_composition([id_report_fp], metadata_fp, clstr_fp, output_fp)

            df = pd.read_table(output_fp, index_col=False)

            self.assertEqual(
                list(df.columns),
                ['CLUSTER', 'ACCESSION', 'SEGMENT', 'GENOTYPE', 'HOST', 'COUNTRY', 'COLLECTION_DATE'],
            )
            self.assertEqual(df['CLUSTER'].tolist(), ['>Cluster 1', '>Cluster 1'])
            self.assertEqual(df['ACCESSION'].tolist(), ['ACC1', 'ACC2'])

    def test_export_cluster_composition_deduplicates_batch_clusters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_fp = os.path.join(tmpdir, 'metadata.csv')
            clstr_fp = os.path.join(tmpdir, 'clusters.clstr')
            id_report_a_fp = os.path.join(tmpdir, 'sample_a_ID_Report.txt')
            id_report_b_fp = os.path.join(tmpdir, 'sample_b_ID_Report.txt')
            output_fp = os.path.join(tmpdir, 'batch_cluster_composition.tsv')

            with open(metadata_fp, 'w', encoding='utf-8') as handle:
                handle.write('ACCESSION;GENOTYPE;SEGMENT;HOST;COUNTRY;COLLECTION_DATE\n')
                handle.write('ACC1;H5N1;4;Gallus gallus;Germany;2024-01-01\n')
                handle.write('ACC2;H5N1;4;Anas platyrhynchos;France;2024-01-02\n')
                handle.write('ACC3;H1N1;6;Homo sapiens;Brazil;2024-01-03\n')

            with open(clstr_fp, 'w', encoding='utf-8') as handle:
                handle.write('>Cluster 1\n')
                handle.write('0\t100nt, >ACC1... *\n')
                handle.write('1\t100nt, >ACC2... at 99%\n')
                handle.write('>Cluster 2\n')
                handle.write('0\t100nt, >ACC3... *\n')

            for fp, cluster in ((id_report_a_fp, '>Cluster 1'), (id_report_b_fp, '>Cluster 1')):
                with open(fp, 'w', encoding='utf-8') as handle:
                    handle.write('SAMPLE_NAME\tSEGMENT\tCLUSTER\tREPRESENTATIVE\t%ID\tASSIGNED_BY\tGENOTYPE\tHOST\tCOUNTRY\n')
                    handle.write(
                        f">sample\t{{'4': 2}}\t{cluster}\tREP1\t99.1\tCD-HIT\t{{'H5N1': 1.0}}\t{{'Gallus gallus': 0.5, 'Anas platyrhynchos': 0.5}}\t{{'Germany': 0.5, 'France': 0.5}}\n"
                    )

            export_cluster_composition([id_report_a_fp, id_report_b_fp], metadata_fp, clstr_fp, output_fp)

            df = pd.read_table(output_fp, index_col=False)

            self.assertEqual(len(df), 2)
            self.assertEqual(df['ACCESSION'].tolist(), ['ACC1', 'ACC2'])
            self.assertEqual(len(df[['CLUSTER', 'ACCESSION']].drop_duplicates()), 2)


if __name__ == '__main__':
    unittest.main()
