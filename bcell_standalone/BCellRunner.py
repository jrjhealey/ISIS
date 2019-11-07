# A runner for Bcell linear epitope prediction
# Because this program is crappy, this has to be run in py2!!!

from antibody_epitope_prediction import AntibodyEpitopePrediction as AEP
from util import print_chart_table
class BCellRunner(object):
    """A class to run, and then hold data from, the various
       B-Cell epitope prediction tools
    """

    def __init__(self, infile, window_size):

        # List methods to enable looping of all method options
        self.methods = ["Chou-Fasman",
                        "Emini",
                        "Karplus-Schulz",
                        "Kolaskar-Tongaonkar",
                        "Parker"]
# These methods require the tcsh interpreter and netsurfp
# and are therefore a pain in the ass
#                   6: "BepiPred-1.0"
#                   7: "BepiPred-2.0"}

# Window size could also be customised, but is set by default in a
# method-specific manner, so left as False
        self.results = {}
        for method in self.methods:
            handle = AEP()
            self.results[method] = handle.predict_antibody_epitope(
            method_name=method, filename=infile,
            swissprot=None, window_size=window_size
            )
            #For debugging etc
            #print(method)
            #print_chart_table(self.results[method])
