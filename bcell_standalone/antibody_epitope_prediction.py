"""
Created on: 03/18/2014
Updated on: 05/19/2017

@author: Dorjee Gyaltsen
@brief: antibody epitope prediction - the main class
"""

import os
import sys
import pickle
from util import *

from collections import OrderedDict
from method_calculation import MethodScores
from collections import namedtuple
SequenceInfo = namedtuple('SequenceInfo', ['sequence_name', 'sequence'])

# adding dependencies to the Python path
project_dir = os.path.dirname(os.path.realpath(__file__))
deps = os.path.join(project_dir, "deps")
sys.path.append(os.path.join(deps, 'iedbtools-utilities'))

from iedbtools_utilities.sequence_io import FASTASequenceInput, OneSequenceInput, WhitespaceSeparatedSequenceInput

class AntibodyEpitopePrediction(object):
    def __init__(self):
        file_path = os.path.join(os.path.dirname(__file__), 'method_scales.pickle')
        pickle_file = open(file_path, 'rb')
        self.scale_dict = pickle.load(pickle_file)

    def predict_antibody_epitope(self, method_name, swissprot, filename, window_size=None):
        if method_name == 'Bepipred-1.0':
            method_name = 'Bepipred'
        elif method_name == 'Bepipred-2.0':
            method_name = 'Bepipred2'

        ms = MethodScores()
        # set values for a particular method
        ms.scale_dict = self.scale_dict[method_name]

        window_size = ms.window if window_size is False else window_size
        center_position = "%d" % round(int(window_size) / 2.0)

        if swissprot:
            sequence_text = get_sequence(swissprot)
        else:
            with open(filename, "r") as infile:
                sequence_text = infile.read()

        if FASTASequenceInput.valid_sequence_text(sequence_text):
            sequence_format = 'fasta'
        elif WhitespaceSeparatedSequenceInput.valid_sequence_text(sequence_text):
            sequence_format = 'whitespace_separated'
        elif OneSequenceInput.valid_sequence_text(sequence_text):
            sequence_format = 'one_sequence'
        else:
            print("Format of sequence not recognized: {}".format(sequence_text))

        # Get the appropriate sequence input class.
        if sequence_format == 'fasta':
            si = FASTASequenceInput(sequence_text)
        elif sequence_format == 'whitespace_separated':
            si = WhitespaceSeparatedSequenceInput(sequence_text)
        elif sequence_format == 'one_sequence':
            si = OneSequenceInput(sequence_text)
        else:
            print("Unexpected sequence_format value: {}".format(sequence_format))

        output_dict = OrderedDict()
        for i, aa_sequence in enumerate(si.as_amino_acid_text()):
            if method_name == 'Emini':
                result = ms.emini_method(method_name, aa_sequence, window_size, center_position)
            elif method_name == 'Karplus-Schulz':
                result = ms.karplusshulz_method(method_name, aa_sequence, window_size, center_position)
            elif method_name == 'Kolaskar-Tongaonkar':
                result = ms.kolaskartongaonkar_method(method_name, aa_sequence, window_size, center_position)
            elif method_name == 'Bepipred':
                result = ms.bepipred_method(method_name, aa_sequence, window_size, center_position)
            elif method_name == 'Bepipred2':
                result = ms.bepipred_method2(method_name, aa_sequence, window_size, center_position)
            else:
                result = ms.classical_method(method_name, aa_sequence, window_size, center_position)
            threshold = round(result[1][0], 3)

            if len(result) == 2:
                epitopes = None
            else:
                epitopes = result[-1]

            # path to the temporary result file
            result_tmpfile = result[0]
            # list of result tuples
            prediction_result = [tuple(row.split(",")) for row in open(result_tmpfile).read().splitlines()]

            # delete temporary result file
            os.unlink(result_tmpfile)

            prediction_input = SequenceInfo(aa_sequence, si.sequence_names[i])

            prediction_output = {
                "epitopes": epitopes,
                "prediction_result": prediction_result,
                "threshold": threshold,
            }

            output_dict.update({
                prediction_input: prediction_output,
            })

        return output_dict
