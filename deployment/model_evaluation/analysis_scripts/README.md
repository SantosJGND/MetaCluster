# Analysis Scripts

Experimental comparison scripts for evaluating model variants on viral
metagenomic data. Each script tests alternative approaches against the
production pipeline and reports side-by-side metrics.

## Setup

```bash
source .venv/bin/activate
export PYTHONPATH=$(pwd)
```

All scripts require prior evaluation pipeline output (study output directory
or cached training data).

---

## Documentation

Detailed descriptions for each script are in the [documentation/](documentation/)
directory:

| Script | Documentation |
|---|---|
| [`compare_sort_strategies.py`](compare_sort_strategies.py) | [compare_sort_strategies.md](documentation/compare_sort_strategies.md) |
| [`composition_model_comparison.py`](composition_model_comparison.py) | [composition_model_comparison.md](documentation/composition_model_comparison.md) |
| [`last_tp_division_prediction_second.py`](last_tp_division_prediction_second.py) | [last_tp_division_prediction.md](documentation/last_tp_division_prediction.md) |
| [`input_composition_prediction.py`](input_composition_prediction.py) | [input_composition_prediction.md](documentation/input_composition_prediction.md) |
| [`debug_metrics_test2.py`](debug_metrics_test2.py) | [debug_metrics_test2.md](documentation/debug_metrics_test2.md) |
| [`../analysis_data_extractor.py`](../analysis_data_extractor.py) | [analysis_data_extractor.md](documentation/analysis_data_extractor.md) |
