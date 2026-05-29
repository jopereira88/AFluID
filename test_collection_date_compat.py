import math
import os
import tempfile
import unittest
from copy import deepcopy
from unittest import mock

from install import _load_cluster_metadata_rows
from main import _normalize_collection_date as normalize_main_collection_date
from main import _unpack_cluster_payload, bclust, redirector
from final_report_utils import _normalize_collection_date as normalize_html_collection_date
from structures import flagdict


class InstallCollectionDateCompatibilityTests(unittest.TestCase):
    def test_load_cluster_metadata_rows_accepts_missing_collection_date(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_fp = os.path.join(tmpdir, "metadata.csv")
            with open(metadata_fp, "w", encoding="utf-8") as handle:
                handle.write("ACCESSION;GENOTYPE;SEGMENT;HOST;COUNTRY\n")
                handle.write("ACC1;H5N1;4;Gallus gallus;Germany\n")

            rows = _load_cluster_metadata_rows(metadata_fp)

        self.assertEqual(rows["ACC1"], ["H5N1", "4", "Gallus gallus", "Germany", ""])


class MainCollectionDateCompatibilityTests(unittest.TestCase):
    def test_unpack_cluster_payload_supports_legacy_shape(self):
        representative, genotypes, segments, hosts, countries, collection_date = _unpack_cluster_payload(
            [{"H5N1": 1}, {"4": 1}, {"Gallus gallus": 1}, {"Germany": 1}, "REP1"]
        )

        self.assertEqual(representative, "REP1")
        self.assertEqual(genotypes, {"H5N1": 1})
        self.assertEqual(segments, {"4": 1})
        self.assertEqual(hosts, {"Gallus gallus": 1})
        self.assertEqual(countries, {"Germany": 1})
        self.assertEqual(collection_date, "Unknown")

    def test_unpack_cluster_payload_supports_new_shape(self):
        representative, _, _, _, _, collection_date = _unpack_cluster_payload(
            [{"H5N1": 1}, {"4": 1}, {"Gallus gallus": 1}, {"Germany": 1}, "2021-01", "REP2"]
        )

        self.assertEqual(representative, "REP2")
        self.assertEqual(collection_date, "2021-01")

    def test_redirector_keeps_rows_with_blank_collection_date(self):
        flags = deepcopy(flagdict)

        with tempfile.TemporaryDirectory() as tmpdir:
            runs_p = tmpdir
            reports_p = tmpdir
            file_tag = "compat"
            mappings = {">Seq_1": ">orig1"}

            with open(os.path.join(runs_p, f"format_{file_tag}.fasta"), "w", encoding="utf-8") as handle:
                handle.write(">Seq_1\n" "ATGCATGCATGCATGCATGC\n")

            with open(os.path.join(reports_p, f"{file_tag}_ID_Report.txt"), "w", encoding="utf-8") as handle:
                handle.write(
                    "SAMPLE_NAME\tSEGMENT\tCLUSTER\tREPRESENTATIVE\t%ID\tASSIGNED_BY\tGENOTYPE\tHOST\tCOUNTRY\tCOLLECTION_DATE\n"
                )
                handle.write(
                    ">orig1\t3\tUNCLUSTERED\tREP1\t98.5\tL-BLAST\tH5N1\tGallus gallus\tGermany\t\n"
                )

            redirector(
                report=f"{file_tag}_ID_Report.txt",
                flags=flags,
                file_tag=file_tag,
                mappings=mappings,
                runs_p=runs_p,
                reports_p=reports_p,
                force_flumut=False,
                force_genin=False,
                force_getref=False,
                mode="consensus",
                single_sample=True,
            )

        self.assertEqual(flags["Sample"]["PA"], [">Seq_1"])
        self.assertEqual(flags["Sample"]["PA_ref"], ["REP1"])

    def test_main_collection_date_normalizer_maps_blank_to_unknown(self):
        self.assertEqual(normalize_main_collection_date("   "), "Unknown")
        self.assertEqual(normalize_main_collection_date(None), "Unknown")

    def test_bclust_accepts_legacy_cluster_metadata_without_collection_date(self):
        flags = deepcopy(flagdict)
        flags["BLAST"]["Sequences unassigned against cluster representatives"] = []

        with tempfile.TemporaryDirectory() as tmpdir:
            samples_p = os.path.join(tmpdir, "samples")
            metadata_p = os.path.join(tmpdir, "metadata")
            blast_p = os.path.join(tmpdir, "blast")
            runs_p = os.path.join(tmpdir, "runs")
            for path in (samples_p, metadata_p, blast_p, runs_p):
                os.makedirs(path, exist_ok=True)

            fasta_name = "to_blast.fasta"
            with open(os.path.join(samples_p, fasta_name), "w", encoding="utf-8") as handle:
                handle.write(">Seq_1\nATGCATGCATGC\n")

            metadata_name = "cluster_desc.txt"
            with open(os.path.join(metadata_p, metadata_name), "w", encoding="utf-8") as handle:
                handle.write(
                    "Cluster\tRepresentative\tGenotypes\tSegments\tHosts\tCountries\n"
                )
                handle.write(
                    ">Cluster 1\tREP1\t{'H5N1': 1}\t{'4': 1}\t{'Gallus gallus': 1}\t{'Germany': 1}\n"
                )

            with mock.patch("main._run_command") as mocked_run, mock.patch(
                "main.best_blast", return_value={"Seq_1": ["REP1", 99.0]}
            ):
                result = bclust(
                    metadata_p=metadata_p,
                    metadata_f=metadata_name,
                    samples_p=samples_p,
                    blast_p=blast_p,
                    blast_db="dummy_db",
                    runs_p=runs_p,
                    file_tag="compat_bclust",
                    flags=flags,
                    fasta=fasta_name,
                )

        mocked_run.assert_called_once()
        self.assertEqual(result["Seq_1"][2].strip(), ">Cluster 1")
        self.assertEqual(result["Seq_1"][7], "Unknown")


class FinalReportCollectionDateCompatibilityTests(unittest.TestCase):
    def test_html_collection_date_normalizer_maps_nan_to_unknown(self):
        self.assertEqual(normalize_html_collection_date(float("nan")), "Unknown")
        self.assertEqual(normalize_html_collection_date("nan"), "Unknown")
        self.assertEqual(normalize_html_collection_date("2022-02-04"), "2022-02-04")


if __name__ == "__main__":
    unittest.main()
