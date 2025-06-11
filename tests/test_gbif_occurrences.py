import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import tempfile
import shutil
import json
from scripts.gbif_occurrences import get_closest_gbif

class TestGBIFOccurrences(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.species = "Test species"
        self.ref_lat = 40.0
        self.ref_lon = -7.5
        self.radius_km = 50
        self.top_n = 1
        self.min_year = 2015

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch("scripts.gbif_occurrences.requests.get")
    def test_get_closest_gbif_success(self, mock_get):
        # Mock GBIF API response with coordinates
        mock_data = {
            "count": 1,
            "results": [
                {
                    "key": 123456,
                    "decimalLatitude": 40.01,
                    "decimalLongitude": -7.49,
                    "eventDate": "2016-05-20",
                    "country": "Portugal",
                    "locality": "Serra da Estrela",
                    "datasetKey": "abcd-1234"
                }
            ]
        }
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = mock_data
        mock_get.return_value = mock_response

        results, no_hit = get_closest_gbif(
            self.species,
            self.ref_lat,
            self.ref_lon,
            self.radius_km,
            self.temp_dir,
            self.top_n,
            self.min_year
        )

        self.assertIsNone(no_hit)
        self.assertIsInstance(results, list)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]["species"], self.species)
        self.assertIn("distance_km", results[0])

    @patch("scripts.gbif_occurrences.requests.get")
    def test_get_closest_gbif_no_occurrence(self, mock_get):
        # Simulate empty result from GBIF API
        mock_data = {"count": 0, "results": []}
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = mock_data
        mock_get.return_value = mock_response

        results, no_hit = get_closest_gbif(
            self.species,
            self.ref_lat,
            self.ref_lon,
            self.radius_km,
            self.temp_dir,
            self.top_n,
            self.min_year
        )

        self.assertIsNone(results)
        self.assertEqual(no_hit, self.species)
