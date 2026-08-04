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
    def register_command(bundle_info, command_info, logger):
        from . import cmd
        cmd.register_command(command_info.name, logger)


bundle_api = _ISISBundle()
