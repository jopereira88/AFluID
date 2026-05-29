import unittest
import importlib
import sys
import types

from metadata_tree_utils import load_reference_dicts


class MetadataTreeUtilsCompatibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geo_dict, cls.taxa_dict = load_reference_dicts(
            "metadata",
            "geo_tree.json",
            "host_taxonomy_tree.json",
        )

    def test_usa_maps_to_legacy_north_america(self):
        self.assertEqual(
            self.geo_dict["United States of America"],
            ("Northern America", "North America"),
        )

    def test_brazil_maps_to_legacy_south_america(self):
        self.assertEqual(
            self.geo_dict["Brazil"],
            ("South America", "South America"),
        )

    def test_viet_nam_keeps_legacy_region_spelling(self):
        self.assertEqual(
            self.geo_dict["Viet Nam"],
            ("South-Eastern Asia", "Asia"),
        )

    def test_subfamily_falls_back_to_family(self):
        self.assertEqual(
            self.taxa_dict["Spatula versicolor"],
            [
                "Spatula versicolor",
                "Spatula",
                "Anatidae",
                "Anatidae",
                "Anseriformes",
                "Anseriformes",
                "Aves",
            ],
        )

    def test_division_falls_back_to_order_then_class(self):
        self.assertEqual(
            self.taxa_dict["Anatidae"],
            [
                "X",
                "X",
                "Anatidae",
                "Anatidae",
                "Anseriformes",
                "Anseriformes",
                "Aves",
            ],
        )
        self.assertEqual(
            self.taxa_dict["Aves"],
            ["X", "X", "X", "X", "Aves", "X", "Aves"],
        )


class FinalReportCompatibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pd = types.ModuleType("pandas")
        pd.isna = lambda value: value is None

        class DataFrame:
            pass

        pd.DataFrame = DataFrame
        sys.modules["pandas"] = pd

        final_report_utils = importlib.import_module("final_report_utils")

        cls.rollup_taxa = staticmethod(final_report_utils.rollup_taxa)
        cls.geo_dict, cls.taxa_dict = load_reference_dicts(
            "metadata",
            "geo_tree.json",
            "host_taxonomy_tree.json",
        )

    def test_perdicinae_rolls_up_via_current_family(self):
        self.assertEqual(self.rollup_taxa("Perdicinae", self.taxa_dict, "Family"), "Phasianidae")
        self.assertEqual(self.rollup_taxa("Perdicinae", self.taxa_dict, "Order"), "Galliformes")
        self.assertEqual(self.rollup_taxa("Perdicinae", self.taxa_dict, "Class"), "Aves")


if __name__ == "__main__":
    unittest.main()
