#!/usr/bin/env python3
"""Export small Hippo/S2 plotting tables, without creating or uploading H5MUs."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cc_perturbseq.compact_plot_inputs import main
if __name__ == "__main__":
    main()
