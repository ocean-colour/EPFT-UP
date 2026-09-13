"""MOANA: reimplementation of the NASA PACE picophytoplankton algorithm.

A separate track from the IOP work (Q&A #3): empirical PC regression from
hyperspectral Rrs (+SST) to Prochlorococcus / Synechococcus / picoeukaryote
cell abundances. Design: ``docs/design/moana_design.md``; background:
``reports/MOANA_Claude_Report.md``.

Public surface:

- :func:`epft_up.moana.algorithm.run_moana` — the retrieval;
- :func:`epft_up.moana.io.load_luts` — the vendored NASA constants;
- :func:`epft_up.moana.pipeline.process_cruise` /
  :func:`epft_up.moana.pipeline.build_training_matrix` — AMT24 Level-2 →
  training matrix;
- :func:`epft_up.moana.train.train_moana` — retraining + basis comparison.
"""

from epft_up.moana.algorithm import run_moana
from epft_up.moana.io import load_luts
from epft_up.moana.pipeline import (DEFAULT_PIPELINE, build_training_matrix,
                                    process_cruise, process_day)
from epft_up.moana.train import compare_loadings, train_moana

__all__ = ['run_moana', 'load_luts', 'process_day', 'process_cruise',
           'build_training_matrix', 'train_moana', 'compare_loadings',
           'DEFAULT_PIPELINE']
