"""File/region/resolution plumbing, including the cold-start twins."""

import cooler
import pytest

from hic_mcp import data
from hic_mcp.data import (
    DataError,
    demo_path,
    list_resolutions,
    open_matrix,
    parse_region_checked,
    resolve_input_path,
)


def test_demo_resolves_and_exists():
    p = resolve_input_path(None)
    assert p.name == "hff_microc_chr17_hg38.mcool"
    assert p.exists()


def test_cold_start_demo_missing_names_the_file(monkeypatch, tmp_path):
    monkeypatch.setenv("HIC_MCP_DATA_DIR", str(tmp_path))
    with pytest.raises(DataError, match="demo dataset is missing"):
        resolve_input_path(None)


def test_user_file_works_without_demo_data(monkeypatch, tmp_path):
    """A user-supplied path must not depend on the bundled data at all."""
    real = demo_path()
    monkeypatch.setenv("HIC_MCP_DATA_DIR", str(tmp_path))
    p = resolve_input_path(str(real))
    assert p == real


def test_missing_file_and_wrong_suffix():
    with pytest.raises(DataError, match="No such file"):
        resolve_input_path("/nonexistent/thing.mcool")
    with pytest.raises(DataError, match="not a .cool"):
        resolve_input_path(__file__)


def test_resolutions_and_bad_resolution():
    p = resolve_input_path(None)
    assert list_resolutions(p) == [10_000, 100_000, 1_000_000]
    with pytest.raises(DataError, match="Available: 10000, 100000, 1000000"):
        open_matrix(p, 25_000, default=10_000)


def test_single_resolution_cool(tiny_unbalanced_cool):
    p = resolve_input_path(tiny_unbalanced_cool)
    assert list_resolutions(p) == [10_000]
    clr = open_matrix(p, None, default=10_000)
    assert clr.binsize == 10_000


def test_region_parse_contract_matches_cooler():
    """Contract-drift guard: any region we accept, cooler's own fetch accepts too."""
    p = resolve_input_path(None)
    clr = open_matrix(p, 100_000, default=100_000)
    for region in ["chr17", "chr17:1,000,000-2,000,000", "chr17:0-100000"]:
        chrom, start, end = parse_region_checked(clr, region)
        m = clr.matrix(balance=False).fetch(region)
        assert m.shape[0] == -(-(end - start) // 100_000)


def test_region_errors_are_agent_readable():
    p = resolve_input_path(None)
    clr = open_matrix(p, 100_000, default=100_000)
    with pytest.raises(DataError, match="Chromosomes here: chr17"):
        parse_region_checked(clr, "chr2:1-1000")
    with pytest.raises(DataError):
        parse_region_checked(clr, "chr17:90,000,000-95,000,000")  # beyond chrom end
    with pytest.raises(DataError, match="zero width"):
        parse_region_checked(clr, "chr17:50,000,000-50,000,000")


def test_data_dir_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("HIC_MCP_DATA_DIR", str(tmp_path))
    assert data.data_dir() == tmp_path


def test_tiny_fixture_is_real_cooler(tiny_unbalanced_cool):
    info = cooler.Cooler(tiny_unbalanced_cool).info
    assert info["nnz"] > 0
