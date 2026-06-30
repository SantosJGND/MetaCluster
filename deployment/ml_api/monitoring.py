import csv
import time
import os
from pathlib import Path
from datetime import datetime

LOG_DIR = Path(__file__).parent.parent / "logs"
LOG_FILE = LOG_DIR / "predictions.csv"


def _ensure_file():
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    if not LOG_FILE.exists():
        with open(LOG_FILE, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["timestamp", "model_version", "prediction", "latency_ms"])


def log_prediction(model_version, prediction, latency_ms):
    _ensure_file()
    with open(LOG_FILE, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([datetime.utcnow().isoformat(), model_version, prediction, round(latency_ms, 2)])


def get_metrics():
    _ensure_file()
    import csv
    count = 0
    total_latency = 0.0
    latest_version = None
    with open(LOG_FILE, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            count += 1
            total_latency += float(row["latency_ms"])
            latest_version = row["model_version"]
    avg_latency = round(total_latency / count, 2) if count > 0 else 0.0
    return {
        "request_count": count,
        "avg_latency_ms": avg_latency,
        "model_version": latest_version,
    }

