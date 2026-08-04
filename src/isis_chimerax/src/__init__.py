"""
ISIS ChimeraX Plugin - B-cell epitope prediction
"""

from chimerax.core.toolshed import BundleAPI

_commands_registered = False


class _ISISBundle(BundleAPI):

    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        """Register commands when bundle initializes."""
        global _commands_registered
        if not _commands_registered:
            from . import cmd
            cmd.register_all_commands(session)
            _commands_registered = True

    @staticmethod
    def start_tool(session, tool_name):
        # Register commands if not already done
        global _commands_registered
        if not _commands_registered:
            from . import cmd
            cmd.register_all_commands(session)
            _commands_registered = True

    @staticmethod
    def finish(session, bundle_info):
        pass


bundle_api = _ISISBundle()
