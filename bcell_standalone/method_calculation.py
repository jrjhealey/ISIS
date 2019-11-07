'''
Created on 03.18.2014
@author: Dorjee Tamang
'''

import os
import tempfile

class MethodScores(object):

    @property
    def scale_dict(self):
        return dict(list(zip(self._dict['residue'], self._dict['score'])))

    @scale_dict.setter
    def scale_dict(self, val):
        self._dict = val

    @property
    def id(self):
        return self._dict['id']

    @property
    def residue(self):
        return self._dict['residue']

    @property
    def score(self):
        return self._dict['score']

    @property
    def window(self):
        return self._dict['default_window']

    @property
    def title(self):
        return self._dict['title']

    @property
    def reference(self):
        return self._dict['reference']

    @property
    def pubmed(self):
        return self._dict['pubmed']

    @property
    def n_scale(self):
        return self._dict['n_scale']

    def classical_method(self, name, sequence, window, center):

        csv_file, stats = self.linear_average(name, self.scale_dict, sequence, int(window), int(center))
        return  csv_file, stats

    def emini_method(self, name, sequence, window, center):

        csv_file, list_values = self.linear_average(name, self.scale_dict, sequence, int(window), int(center))

        stats = self.compute_statistics(list_values)
        threshold = 1
        min_pred_size = 6
        return csv_file, stats, self.predict_epitope(list_values, min_pred_size, threshold)

    def karplusshulz_method(self, name, sequence, window, center):

        AA = ["K", "S", "G", "P", "D", "E", "Q", "T", "N", "R", "A", "L", "H", "V", "Y", "I", "F", "C", "W", "M"]
        BNORM0 = [1.093,1.169,1.142,1.055,1.033,1.094,1.165,1.073,1.117,1.038,1.041,0.967,0.982,0.982,0.961,1.002,0.930,0.960,0.925,0.947]
        BNORM1= [1.082,1.048,1.042,1.085,1.089,1.036,1.028,1.051,1.006,1.028,0.946,0.961,0.952,0.927,0.930,0.892,0.912,0.878,0.917,0.862]
        BNORM2 = [1.057,0.923,0.923,0.932,0.932,0.933,0.885,0.934,0.930,0.901,0.892,0.921,0.894,0.913,0.837,0.872,0.914,0.925,0.803,0.804]
        WT = [0.25,0.50,0.75,1.00,0.75,0.50,0.25]

        nRes = len(sequence)
        NAYB = []
        # Original
        for i in range(1, len(sequence)-1):
        # My edit:
        #for i in range(0, len(sequence)):
            val = 0
            # Original
            i1 = AA.index(sequence[i-1:i])
            i2 = AA.index(sequence[i+1:i+2])
            # My edit
        #    i1 = AA.index(sequence[i-1:i])
        #    i2 = AA.index(sequence[i+1:i+2])

            if i1 >= 10 or i2 >= 10:
                val = 1

            if i1 >= 10 and i2 >= 10:
                val = 2

            NAYB.append(val)

        # TODO: check if this is the right way
        NAYB.insert(0, 0)
        NAYB.append(0)

        list_values = []
        for i in range(4, len(sequence)-3):
            sum = 0
            peptide = ""
            for j in range(0, int(window)):
                res = sequence[i-4+j:i-4+j+1]
                peptide += res
                index = AA.index(res)

                if NAYB[i-4+j] == 0:
                    sum += BNORM0[index] * WT[j] / 4.0
                if NAYB[i-4+j] == 1:
                    sum += BNORM1[index] * WT[j] / 4.0
                if NAYB[i-4+j] == 2:
                    sum += BNORM2[index] * WT[j] / 4.0

            position = i
            residue = sequence[position-1:position]
            startPos = i - 4 + 1
            endPos = startPos + int(window) - 1
            list_values.append([position, residue, startPos, endPos, peptide, sum])

        # save the values into a csv file
        csv_file = self.create_csv(list_values)

        stats = self.compute_statistics(list_values)
        threshold = 1
        return  csv_file, stats

    def kolaskartongaonkar_method(self, name, sequence, window, center):

        cutoff = 1
        csv_file, list_values = self.linear_average(name, self.scale_dict, sequence, int(window), int(center))
        stats = self.compute_statistics(list_values)
        #TODO: calculate the threshold
        threshold = 1
        min_pred_size = 6
        return csv_file, stats, self.predict_epitope(list_values, min_pred_size, threshold)


    def bepipred_method(self, name, sequence, window, center):
        threshold = 0
        # create a temporary directory (if doesn't exist)
        tmpfile = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
        tmpfile.write(">sequence 1\n%s\n" %sequence)
        tmpfile.close()
        return self.predict_bepipred(tmpfile, sequence, threshold)

    def bepipred_method2(self, name, sequence, window, center, threshold = 0):
        #threshold = 0
        # create a temporary directory (if doesn't exist)
        tmpfile = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
        tmpfile.write(">sequence 1\n%s\n" %sequence)
        tmpfile.close()
        return self.predict_bepipred2(tmpfile, sequence, threshold)

    def linear_average(self, name, scale, sequence, window, center):
        import math

        nRes = len(sequence)
        list_values= []
        for i in range(nRes+1 - window):
            product = 1
            total_scale = 0
            peptide = ""

            for j in range(window):
                res = sequence[i + j : i + j + 1]
                peptide += res
                if name == "Emini": product *= scale[res]
                else: total_scale += scale[res]

            if name == "Emini": product *= math.pow(0.37, -6)
            else: average = total_scale / window

            position = i + center
            residue = sequence[position - 1: position]
            start_position = i + 1
            end_position = start_position + window - 1
            opt = product if name == "Emini" else average
            list_values.append([position, residue, start_position, end_position, peptide, opt])

        stats = self.compute_statistics(list_values)

        if name == "Emini":
            list_values = self.normalize(list_values)

        # save the values into a csv file
        tmpfile_name = self.create_csv(list_values)

        if name == "Emini" or name == "Kolaskar-Tongaonkar":
            return tmpfile_name, list_values
        else:
            return tmpfile_name, stats

    def compute_statistics(self, result, opt=None, value=None):
        if not value:
            position = []
            value = []

            for res in result:
                position.append(res[0])
                if not opt:
                    value.append(res[5])
                else: value.append(res[2])

        average = sum(map(float, value))/float(len(value))
        min_value = min(value)
        max_value = max(value)
        min_index = value.index(min(value))
        max_index = value.index(max(value))

        return average, min_value, max_value, min_index, max_index

    def create_csv(self, list_values):
        import csv

        # bcell_tmpdir = './output'

        # create a temporary file inside the tmp/ directory
        tmpfile = tempfile.NamedTemporaryFile(prefix="csv_", delete=False)

        with open(tmpfile.name, 'wb') as result:
            writer = csv.writer(result)

            if self.id != 'Bepipred':
                writer.writerow(('Position','Residue','Start','End','Peptide','Score'))
            else: writer.writerow(('Position','Residue','Score','Assignment'))

            for values in list_values:
                if isinstance(values[-1], float):
                    values[-1] = round(values[-1], 3)
                writer.writerow(values)

        tmpfile.close()
        return tmpfile.name

    def normalize(self, list_values):
        stats = self.compute_statistics(list_values)
        average = stats[0]
        list_normalize = []
        for idx, lis in enumerate(list_values):
            lis[5] = lis[5]/average
            list_normalize.append(lis)
        return list_normalize

    def predict_epitope(self, result, min, threshold):
        size = 0
        start_indx = 0
        seq = ""
        average = 0
        bFirst = False
        pred_peptide = []
        for idx, ele in enumerate(result):
            if ele[5] >= threshold:
                size += 1
                seq += ele[1]
                average += ele[5]
                if not bFirst:
                    start_indx = idx
                    bFirst = True
            else:
                if size >= min:
                    start_pos = ele[0]-len(seq)
                    end_pos = ele[0]-1
                    pred_peptide.append((start_pos, end_pos, seq, average/len(seq)))

                size = 0
                start_indx = idx
                seq = ""
                average = 0
                bFirst = False

        return pred_peptide

    def predict_bepipred(self, infile, sequence, threshold):
        import subprocess as sub

        bepipred_executable_dir = "bepipred-1.0b"
        path_method = os.path.join(os.path.dirname(__file__), bepipred_executable_dir)
        path_executable = os.path.join(path_method, 'bepipred')

        proc = sub.Popen([path_executable, '-t', str(threshold) , infile.name], stdout=sub.PIPE, stderr=sub.PIPE)
        output, error = proc.communicate()
        if error:
            msg = "error occurred!"

        lines = output.splitlines(True)
        new_list = [line.strip() for line in lines if not line.startswith("#")]

        pos = []
        score = []
        assignment = []

        for ln in new_list:
            pos.append(ln.split()[3])
            score.append(ln.split()[5])
            assignment.append(ln.split()[-1])

        pred_peptides = self.get_predicted_peptides(assignment, sequence)

        list_values = list(zip(pos, list(sequence), score, assignment))
        list_values = list(list(x) for x in list_values)

        stats = self.compute_statistics(list_values, True)

        tmpfile_name = self.create_csv(list_values)

        # remove temp file from /tmp directory
        os.remove(infile.name)
        return tmpfile_name, stats, self.get_predicted_peptides(assignment, sequence)

    def predict_bepipred2(self, infile, sequence, threshold):
        import bepipred2 as bp2 #imports BepiPred-2.0
        bp2.utils.RF_MODEL = bp2.utils.init_rf() #Unloads the random forest predictor

        seq = 'cdafvgtwKLVssenfddymkevgvgfatrkvagMAKpnmiisvngdlvtirsesTfkn' #Amino acid sequence
        seq = sequence
        id  = 'Example1' #Name/identifier of sequence
        AG = bp2.Antigen(id, seq) #Sets up antigen class object
        AG.pred_netsurfp() #Predicts surface and secondary structure of antigen
        AG.get_features() #Sets up the feature space for predicting
        AG.predict() #Predicts the epitopes of the antigen
        scores = AG.predicted
        pos = []
        score = []
        assignment = []

        for i in range(len(scores)):
            pos.append(str(i))
            score.append('%.3f' % scores[i])
            if scores[i] >= threshold:
                assignment.append('E')
            else:
                assignment.append('.')

        list_values = list(zip(pos, list(sequence), score, assignment))
        list_values = list(list(x) for x in list_values)

        stats = self.compute_statistics(list_values, True)

        tmpfile_name = self.create_csv(list_values)

        # remove temp file from /tmp directory
        os.remove(infile.name)
        return tmpfile_name, stats, self.get_predicted_peptides(assignment, sequence)

    def get_predicted_peptides(self, assignment, sequence):
        peptides = []
        epiNo = 0
        epiStart = 0
        epiEnd = 0

        for i, ass in enumerate(assignment):
            peptide = ""
            if i == 0 and assignment[i] == "E":
                peptide += sequence[i]
                epiStart = i
            else:
                if i > 0:
                    if assignment[i] == "E":
                        if assignment[i-1] == ".":
                            peptide += sequence[i]
                            epiStart = i
                        else:
                            peptide += sequence[i]
                    else:
                        if assignment[i] == "." and assignment[i-1] == "E":
                            epiEnd = i-1
                            epiNo += 1

                            peptides.append((epiNo, epiStart+1, epiEnd+1, sequence[epiStart:epiEnd+1], (epiEnd - epiStart)+1))
        return peptides
