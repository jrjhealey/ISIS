"""
ISIS ChimeraX Plugin - B-cell epitope prediction
"""

from chimerax.core.toolshed import BundleAPI


class _ISISBundle(BundleAPI):

    api_version = 1

    @staticmethod
    def start_tool(session, tool_name):
        pass

    @staticmethod
    def run_provider(session, name, mgr, **kw):
        """Called when command provider is needed."""
        if mgr == session.toolshed.triggers.command_manager():
            from .cmd import register_command
            register_command(session, name)

    @staticmethod
    def get_class(class_name):
        return None

    @staticmethod
    def initialize(session, bundle_info):
        """Register commands when bundle loads."""
        from .cmd import register_all_commands
        register_all_commands(session)


bundle_api = _ISISBundle()
