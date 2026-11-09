import os
import random
import re
import tempfile
import unittest
from decimal import Decimal

import pandas as pd

from final_report_utils import (
    _normalize_distribution_for_export,
    _serialize_distribution,
    clean_country,
    clean_genotype,
    export_cluster_composition,
    export_id_report_rollup,
    finalize_id_report,
)
from metadata_tree_utils import load_reference_dicts

ID_REPORT_HEADER = 'SAMPLE_NAME\tSEGMENT\tCLUSTER\tREPRESENTATIVE\t%ID\tASSIGNED_BY\tGENOTYPE\tHOST\tCOUNTRY\tCOLLECTION_DATE\n'
ID_REPORT_CLUSTER_ROW = (
    ">sample1\t{'4': 4003}\t>Cluster 10\tREP1\t99.1\tCD-HIT\t{'H1N1': 0.999333, 'H1N2': 0.000667}\t"
    "{'Homo sapiens': 0.96, 'Gallus gallus': 0.001, 'Sus scrofa': 0.038}\t"
    "{'USA': 0.381, 'France': 0.002, 'Reunion': 0.004, 'Brazil': 0.616}\t2009\n"
)
ID_REPORT_LBLAST_ROW = '>sample2\t2\tUNCLUSTERED\tCY000002.1\t99\tL-BLAST\tH1N1\tHomo sapiens\tBrazil\t2024-01\n'


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

    def assertValidDistribution(self, serialized):
        """The written text has at most 2 dp and its values total exactly 1."""
        values = re.findall(r"': ([0-9.eE+-]+)", serialized)
        self.assertTrue(values, serialized)
        for value in values:
            self.assertRegex(value, r'^\d+(\.\d{1,2})?$', serialized)
        self.assertEqual(sum(Decimal(value) for value in values), Decimal(1), serialized)

    def test_distributions_always_have_two_decimals_and_sum_to_one(self):
        rng = random.Random(1988)
        for _ in range(2000):
            size = rng.randint(1, 80)
            counter = {f'k{i}': rng.choice([rng.random(), rng.randint(1, 5000), rng.random() / 1000]) for i in range(size)}
            self.assertValidDistribution(str(_normalize_distribution_for_export(counter)))
            self.assertValidDistribution(_serialize_distribution(counter))

    def test_zero_means_below_half_a_percent(self):
        rng = random.Random(2026)
        for _ in range(2000):
            size = rng.randint(1, 60)
            counter = {f'k{i}': rng.choice([rng.randint(1, 5), rng.randint(1, 500)]) for i in range(size)}
            total = sum(counter.values())
            rounded = _normalize_distribution_for_export(counter)
            self.assertValidDistribution(str(rounded))
            for key, value in rounded.items():
                if value == 0:
                    self.assertLess(counter[key] / total, 0.005, counter)

    def test_shares_above_half_a_percent_are_not_shown_as_zero(self):
        # Larger remainders used to take the spare unit from Colombia (Cluster 21).
        self.assertGreater(_normalize_distribution_for_export({'a': 361, 'b': 122, 'Colombia': 3})['Colombia'], 0)
        # More such entries than spare units: the largest entry gives one up.
        rounded = _normalize_distribution_for_export({'a': 0.006, 'b': 0.006, 'c': 0.008, 'd': 0.98})
        self.assertEqual(rounded, {'d': 0.97, 'c': 0.01, 'a': 0.01, 'b': 0.01})

    def test_distribution_ignores_values_that_cannot_be_proportions(self):
        counter = {'a': 1, 'b': float('nan'), 'c': float('inf'), 'd': -1, 'e': 0}
        self.assertEqual(_normalize_distribution_for_export(counter), {'a': 1.0, 'e': 0.0})

    def test_sum_just_over_one_is_not_treated_as_a_percentage(self):
        cell = "{'USA': 0.504, 'Brazil': 0.501}"
        self.assertEqual(clean_genotype(cell), {'USA': 0.5, 'Brazil': 0.5})
        self.assertEqual(clean_country(cell, self.geo_dict), {'United States of America': 0.5, 'Brazil': 0.5})

    def test_percentages_are_still_normalized(self):
        self.assertEqual(clean_genotype("{'H5N1': '75%', 'H1N1': 25}"), {'H5N1': 0.75, 'H1N1': 0.25})

    def test_rollup_keeps_share_of_entries_that_round_to_zero_individually(self):
        countries = ['Kenya', 'Uganda', 'Ethiopia', 'Mozambique', 'Zambia', 'Zimbabwe', 'Lesotho', 'Botswana']
        cell = {country: 0.004 for country in countries}
        cell['USA'] = 1 - 0.004 * len(countries)
        self.assertEqual(clean_country(cell, self.geo_dict, mode='continent')['Africa'], 0.03)

    def test_rollup_export_rounds_source_columns_and_keeps_scalars(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            id_report_fp = os.path.join(tmpdir, 'sample_ID_Report.txt')
            output_fp = os.path.join(tmpdir, 'sample_ID_Report_rollup.tsv')
            with open(id_report_fp, 'w', encoding='utf-8') as handle:
                handle.write(ID_REPORT_HEADER)
                handle.write(ID_REPORT_CLUSTER_ROW)
                handle.write(ID_REPORT_LBLAST_ROW)

            export_id_report_rollup(id_report_fp, output_fp, self.geo_dict, self.taxa_dict)
            df = pd.read_table(output_fp, index_col=False)

            for col in ['GENOTYPE', 'HOST', 'COUNTRY', 'HOST_FAMILY', 'HOST_CLASS', 'COUNTRY_REGION', 'COUNTRY_CONTINENT']:
                self.assertValidDistribution(df.loc[0, col])
            self.assertEqual(df.loc[1, 'GENOTYPE'], 'H1N1')
            self.assertEqual(df.loc[1, 'COUNTRY'], 'Brazil')

    def test_finalize_id_report_rounds_only_distributions_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            id_report_fp = os.path.join(tmpdir, 'sample_ID_Report.txt')
            with open(id_report_fp, 'w', encoding='utf-8') as handle:
                handle.write(ID_REPORT_HEADER)
                handle.write(ID_REPORT_CLUSTER_ROW)
                handle.write(ID_REPORT_LBLAST_ROW)
                handle.write('>sample3\tNA\tUnassigned\tNA\tNA\tNA\tNA\tNA\tNA\tUnknown\n')

            finalize_id_report(id_report_fp)
            with open(id_report_fp, encoding='utf-8') as handle:
                first_pass = handle.read()
            finalize_id_report(id_report_fp)
            with open(id_report_fp, encoding='utf-8') as handle:
                self.assertEqual(handle.read(), first_pass)

            lines = first_pass.splitlines()
            self.assertEqual(lines[0] + '\n', ID_REPORT_HEADER)
            cells = lines[1].split('\t')
            self.assertEqual(cells[1], "{'4': 4003}")
            for index in (6, 7, 8):
                self.assertValidDistribution(cells[index])
            self.assertEqual(lines[2] + '\n', ID_REPORT_LBLAST_ROW)
            self.assertEqual(lines[3], '>sample3\tNA\tUnassigned\tNA\tNA\tNA\tNA\tNA\tNA\tUnknown')

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
