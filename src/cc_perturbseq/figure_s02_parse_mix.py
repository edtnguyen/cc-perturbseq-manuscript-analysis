"""Figure S2 mixed hiPSC–K562 Parse-v3 data and panel producers."""
from . import figure_s02_parse_mix_data as _data
from .figure_s02_parse_mix_preflight import preflight_parse_mix

# ``read_expression_and_qc`` and ``prepare_parse_mix_benchmark`` were defined in
# the data module, so their global lookup resolves through that module's
# namespace. Replace the legacy preflight there before re-exporting its public
# functions. The legacy implementation incorrectly required ``DGE_filtered``
# itself to be a file even though it is the source notebook's input directory.
_data.preflight_parse_mix = preflight_parse_mix

from .figure_s02_parse_mix_data import *  # noqa: E402,F401,F403
from .figure_s02_parse_mix_data import __all__ as _data_all  # noqa: E402
from .figure_s02_parse_mix_plots import *  # noqa: E402,F401,F403
from .figure_s02_parse_mix_plots import __all__ as _plot_all  # noqa: E402

__all__ = [*_data_all, *_plot_all]
