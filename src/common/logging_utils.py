"""PHI-safe logging helpers.

Rules
-----
- Never log record contents. Log counts, paths, tool versions, exit codes.
- Sample IDs are PHI. Do not print them, even in DEBUG.
- Anything logged from a worker goes to a per-chunk file under $OUTPUT_ROOT/logs/,
  which is gitignored and stays on Minerva.
"""
from __future__ import annotations

import logging
import os
import re
from logging import Logger
from pathlib import Path

# Redact anything that looks like an MSM participant ID.
_ID_PATTERNS = [
    re.compile(r"SINAI_\d+_[A-Z]{2}\d+", re.I),
    re.compile(r"SINAI-Million_\d+_\d+", re.I),
    re.compile(r"\bMRN[0-9]{5,}\b", re.I),
]


class PhiRedactingFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        msg = record.getMessage()
        for pat in _ID_PATTERNS:
            msg = pat.sub("[REDACTED_ID]", msg)
        record.msg = msg
        record.args = ()
        return True


def get_logger(name: str, log_dir: str | os.PathLike | None = None) -> Logger:
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    sh.addFilter(PhiRedactingFilter())
    logger.addHandler(sh)
    if log_dir is not None:
        Path(log_dir).mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(Path(log_dir) / f"{name}.log")
        fh.setFormatter(fmt)
        fh.addFilter(PhiRedactingFilter())
        logger.addHandler(fh)
    return logger
