#!/usr/bin/env python3
"""Backward-compatible entry point for the original script name."""

from map_kmers_ac import main


if __name__ == "__main__":
    raise SystemExit(main())
