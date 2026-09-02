"""Shared fixtures: offline guard and a real (tiny) unbalanced cooler file."""

import socket

import cooler
import numpy as np
import pandas as pd
import pytest


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Every test runs with sockets disabled - the bundled demo needs zero network."""

    def _blocked(*args, **kwargs):
        raise RuntimeError("network access attempted during an offline test")

    monkeypatch.setattr(socket.socket, "connect", _blocked)
    yield


@pytest.fixture(scope="session")
def tiny_unbalanced_cool(tmp_path_factory) -> str:
    """A real single-resolution .cool with no weight column (built with cooler itself)."""
    n, binsize = 120, 10_000
    bins = pd.DataFrame(
        {
            "chrom": "chrT",
            "start": np.arange(n) * binsize,
            "end": (np.arange(n) + 1) * binsize,
        }
    )
    rng = np.random.default_rng(7)
    i, j = np.triu_indices(n)
    counts = rng.poisson(120.0 / (1 + (j - i)) ** 0.9).astype(int)
    keep = counts > 0
    pixels = pd.DataFrame({"bin1_id": i[keep], "bin2_id": j[keep], "count": counts[keep]})
    path = str(tmp_path_factory.mktemp("fixtures") / "tiny_unbalanced.cool")
    cooler.create_cooler(path, bins, pixels)
    return path


@pytest.fixture(scope="session")
def two_chromosome_cool(tmp_path_factory) -> str:
    """A real balanced cooler with TWO chromosomes.

    The bundled demo has one chromosome, so single-chromosome assumptions in the code
    are invisible to every test that uses it - this fixture is the twin that catches
    them, on the bring-your-own-file path the README advertises.
    """
    n1, n2, binsize = 60, 40, 10_000
    bins = pd.DataFrame(
        {
            "chrom": ["cA"] * n1 + ["cB"] * n2,
            "start": list(np.arange(n1) * binsize) + list(np.arange(n2) * binsize),
            "end": list((np.arange(n1) + 1) * binsize) + list((np.arange(n2) + 1) * binsize),
        }
    )
    rng = np.random.default_rng(3)
    i, j = np.triu_indices(len(bins))
    counts = rng.poisson(150.0 / (1 + abs(j - i)) ** 0.8).astype(int)
    keep = counts > 0
    pixels = pd.DataFrame({"bin1_id": i[keep], "bin2_id": j[keep], "count": counts[keep]})
    path = str(tmp_path_factory.mktemp("fixtures") / "two_chrom.cool")
    cooler.create_cooler(path, bins, pixels)
    cooler.balance_cooler(cooler.Cooler(path), store=True, cis_only=True)
    return path
