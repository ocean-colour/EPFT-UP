"""Shared pytest skip guards and hang guard for the EPFT-UP test suite.

Ported from IOPtics ``ioptics/tests/conftest.py`` (branch ``moana`` @
``3aa3b6e``, 2026-09-13), keeping only the idiom and the markers the MOANA
track uses. Two tiers:

- **Tier 1 — data-independent** tests run everywhere (laptop, CI runner) on
  tiny synthetic fixtures and the vendored LUTs in ``epft_up/data/``.
- **Tier 2 — data-dependent** tests skip automatically when the ``$OS_COLOR``
  data tree (or a specific dataset) is unavailable, so the suite is green and
  fast where no data is mounted.

Tier-2 tests opt in with the markers exported here::

    from epft_up.tests.conftest import needs_amt24

    @needs_amt24
    def test_process_day_real():
        ...

The availability probes never raise, so test *collection* cannot fail on a
missing data tree — the guarded tests simply skip.

**Hang guard.** Every test runs under a wall-clock ceiling (see
:func:`_guard_against_hangs`) so a wedged test — a blocked data load, a
download that never returns — fails fast with a traceback instead of hanging
the whole run. It is dependency-free (Unix ``SIGALRM``) and defers to
``pytest-timeout`` when that plugin is installed.
"""

import os
import faulthandler
import signal
import threading

import pytest


# --------------------------------------------------------------------
# Hang guard — a per-test wall-clock ceiling so nothing hangs the run
# --------------------------------------------------------------------
def _pytest_timeout_installed():
    """True if the ``pytest-timeout`` plugin is importable (then defer to it)."""
    try:
        import pytest_timeout       # noqa: F401
        return True
    except Exception:
        return False


def _test_timeout_seconds():
    """Per-test ceiling in seconds; override with ``$EPFT_UP_TEST_TIMEOUT``.

    Defaults to 120 s — far above any real test here (the whole suite runs in
    a few seconds) but low enough that a genuine hang is caught quickly.
    Set to ``0`` (or a negative value) to disable the guard.
    """
    raw = os.getenv('EPFT_UP_TEST_TIMEOUT', '120')
    try:
        return float(raw)
    except (TypeError, ValueError):
        return 120.0


def pytest_configure(config):
    """Register the ``timeout`` marker (also used by ``pytest-timeout``)."""
    config.addinivalue_line(
        'markers',
        'timeout(seconds): override the per-test hang-guard ceiling '
        '(0 disables it for that test).')


@pytest.fixture(autouse=True)
def _guard_against_hangs(request):
    """Fail a test that exceeds the wall-clock ceiling instead of hanging.

    Uses ``SIGALRM``/``setitimer`` (Unix, main thread only). When
    ``pytest-timeout`` is installed we defer to it — it handles some C-level
    hangs this fixture cannot. A hang dumps all thread tracebacks (via
    ``faulthandler``) before failing, so the culprit is easy to spot.
    """
    # Defer to pytest-timeout when present; skip where SIGALRM is unavailable
    # (non-Unix) or unusable (worker thread — signals only fire on the main one).
    if (_pytest_timeout_installed()
            or not hasattr(signal, 'SIGALRM')
            or threading.current_thread() is not threading.main_thread()):
        yield
        return

    timeout = _test_timeout_seconds()
    marker = request.node.get_closest_marker('timeout')
    if marker is not None and marker.args:
        timeout = float(marker.args[0])
    if timeout <= 0:
        yield
        return

    def _on_alarm(signum, frame):
        faulthandler.dump_traceback()
        raise TimeoutError(
            f'test exceeded the {timeout:g}s hang-guard ceiling '
            f'($EPFT_UP_TEST_TIMEOUT); treated as a hang')

    previous = signal.signal(signal.SIGALRM, _on_alarm)
    signal.setitimer(signal.ITIMER_REAL, timeout)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


# --------------------------------------------------------------------
# Tier-2 availability probes and the markers built on them
# --------------------------------------------------------------------
def _os_color_available():
    """True if the ``$OS_COLOR`` data tree is set and present on disk."""
    root = os.getenv('OS_COLOR')
    return root is not None and os.path.isdir(root)


def _amt24_available():
    """True if the AMT24 HyperSAS Level-2 tree is mounted under ``$OS_COLOR``.

    Both layouts are accepted — the tree moved once already
    (``$OS_COLOR/AMT24/`` → ``$OS_COLOR/AMT/AMT24/``, 2026-08-16).
    """
    root = os.getenv('OS_COLOR')
    if root is None:
        return False
    return any(os.path.isdir(os.path.join(root, *parts, 'Radiometry', 'level2'))
               for parts in (('AMT', 'AMT24'), ('AMT24',)))


needs_data = pytest.mark.skipif(
    not _os_color_available(), reason='requires the $OS_COLOR data tree')

needs_amt24 = pytest.mark.skipif(
    not _amt24_available(),
    reason='requires the AMT24 tree under $OS_COLOR (MOANA track)')

#: Earthdata-credentialled tests (MOANA validation target (iii)).
needs_netrc = pytest.mark.skipif(
    not os.path.isfile(os.path.expanduser('~/.netrc')),
    reason='requires Earthdata credentials in ~/.netrc')
