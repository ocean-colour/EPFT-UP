"""Smoke test: the package imports and declares a version."""

import epft_up


def test_version():
    assert isinstance(epft_up.__version__, str)
    assert epft_up.__version__
