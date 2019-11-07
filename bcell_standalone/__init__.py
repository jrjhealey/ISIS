import tempfile
from util import *
from antibody_epitope_prediction import AntibodyEpitopePrediction


def predict(method_name=None, swissprot=None, sequences=None, window_size=None):

    filename = None
    if sequences:
        # check if input is a list
        assert type(sequences) is list, "Input file must be a list of sequence(s)."

        # write a temporary file from a sequence_list items
        tmpfile = tempfile.NamedTemporaryFile()
        for sequence in sequences:
            tmpfile.write("{}\n".format(sequence))
        tmpfile.seek(0)
        filename = tmpfile.name

    aep = AntibodyEpitopePrediction()
    result = aep.predict_antibody_epitope(method_name=method_name, swissprot=swissprot, filename=filename,
                                          window_size=window_size)
    if filename:
        tmpfile.close()

    if not result:
        msg = 'Error calling bcell standalone:\n{}'.format(method_name)
        raise Exception(msg)

    return result
