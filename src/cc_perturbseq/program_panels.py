from .program_io import *
from .program_figure2 import *
from .program_annotation import *
from .program_s4 import *
from .screen_panels import *
from .single_cell_panels import *
from .focused_panels import *
from .qc_panels import *
from .benchmark_panels import *
from .multiome_panels import *
from .panel_adapters import *

from . import multiome_panels as _multiome_panels
from .motif_sources import canonical_enhancer_enrichment_path


def standardize_enrichment(path, *args, **kwargs):
    """Load motif enrichment only from the manuscript-authorized table."""
    canonical = canonical_enhancer_enrichment_path(path)
    return _multiome_panels.standardize_enrichment(
        canonical, *args, **kwargs
    )
