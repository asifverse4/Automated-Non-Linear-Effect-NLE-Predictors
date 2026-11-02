"""
Entrypoint module for running nle_predictor directly via `python -m nle_predictor`.
"""

import sys

from nle_predictor.cli import main

if __name__ == "__main__":
    main(sys.argv[1:])

