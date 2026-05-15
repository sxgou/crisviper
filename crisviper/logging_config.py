"""crisviper/logging_config.py — Logging configuration for crisviper"""

import logging
import sys


# Module-level logger lookup
def get_logger(name: str) -> logging.Logger:
    """Get a logger for the given module name."""
    return logging.getLogger(f"crisviper.{name}")


def setup_logging(verbose: bool = False, quiet: bool = False) -> None:
    """
    Configure the crisviper logger.

    Args:
        verbose: Enable DEBUG level output
        quiet: Suppress INFO and below, show only WARNING/ERROR
    """
    root = logging.getLogger("crisviper")
    root.handlers.clear()

    handler = logging.StreamHandler(sys.stderr)

    if quiet:
        level = logging.WARNING
        formatter = logging.Formatter("%(levelname)s: %(message)s")
    elif verbose:
        level = logging.DEBUG
        formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%H:%M:%S"
        )
    else:
        level = logging.INFO
        formatter = logging.Formatter("%(message)s")  # Plain like print() for info

    handler.setFormatter(formatter)
    root.setLevel(level)
    handler.setLevel(level)
    root.addHandler(handler)

    # Suppress noisy library loggers
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
