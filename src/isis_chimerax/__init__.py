"""
ISIS ChimeraX Plugin

Provides commands for B-cell epitope prediction and visualization
directly within ChimeraX.

Commands:
    isis predict <model> [method <name>] [window <int>] [threshold <float>]
    isis color <model> [method <name>] [palette <colors>]
    isis clear <model>
    isis list
"""

from chimerax.core.toolshed import BundleAPI


class _ISISBundle(BundleAPI):

    api_version = 1

    @staticmethod
    def register_command(bi, ci, logger):
        from . import cmd
        cmd.register_commands(ci, logger)

    @staticmethod
    def get_class(class_name):
        if class_name == "ISISResultsModel":
            from .results import ISISResultsModel
            return ISISResultsModel
        return None


bundle_api = _ISISBundle()
