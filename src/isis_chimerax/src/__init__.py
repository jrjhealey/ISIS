"""
ISIS ChimeraX Plugin - B-cell epitope prediction
"""

from chimerax.core.toolshed import BundleAPI


class _ISISBundle(BundleAPI):

    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        """Register commands when bundle initializes."""
        from . import cmd
        cmd.register_all_commands(session)

    @staticmethod
    def finish(session, bundle_info):
        pass


bundle_api = _ISISBundle()


def get_api():
    return bundle_api
