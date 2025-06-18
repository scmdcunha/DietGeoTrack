from unittest.mock import patch, MagicMock
from scripts.fetch_taxonomy_ncbi import fetch_taxonomy_wrapper

@patch("scripts.fetch_taxonomy_ncbi.get_organism_name", return_value="Eratigena saeva")
@patch("scripts.fetch_taxonomy_ncbi.ncbi.get_taxid_translator")
@patch("scripts.fetch_taxonomy_ncbi.ncbi.get_rank")
@patch("scripts.fetch_taxonomy_ncbi.ncbi.get_lineage")
@patch("scripts.fetch_taxonomy_ncbi.ncbi.get_name_translator")
def test_fetch_taxonomy_wrapper(
    mock_get_name_translator,
    mock_get_lineage,
    mock_get_rank,
    mock_get_translator,
    mock_get_organism_name
):
    # Mock taxonomic info
    mock_get_name_translator.return_value = {"Eratigena saeva": [12345]}
    mock_get_lineage.return_value = [2, 3, 4]
    mock_get_rank.return_value = {2: "order", 3: "family", 4: "genus"}
    mock_get_translator.return_value = {2: "Araneae", 3: "Agelenidae", 4: "Eratigena"}

    result = fetch_taxonomy_wrapper(("hap1234", "ABC123.1", "test@example.com", None))

    assert result is not None, "fetch_taxonomy_wrapper returned None"
    assert result["Query ID"] == "hap1234"
    assert result["Accession ID"] == "ABC123.1"
    assert result["Order"] == "Araneae"
    assert result["Family"] == "Agelenidae"
    assert result["Genus"] == "Eratigena"
    assert result["Species"] == "Eratigena saeva"
