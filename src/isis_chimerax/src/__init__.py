"""
ISIS ChimeraX Plugin - B-cell epitope prediction
"""

from chimerax.core.toolshed import BundleAPI

class _ISISBundle(BundleAPI):

    api_version = 1

    @staticmethod
    def run_provider(session, name, mgr):
        """Register commands when provider is invoked."""
        from .cmd import register_command
        register_command(session, name)

bundle_api = _ISISBundle()
