#!/usr/bin/env python3
"""
Generate method_scales.pickle for antibody epitope prediction.
Run this script once to create the pickle file.
"""

import pickle
import os

AMINO_ACIDS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

SCALES = {
    'Emini': {
        'id': 'Emini',
        'title': 'Emini Surface Accessibility',
        'residue': AMINO_ACIDS,
        'score': [0.815, 1.475, 1.296, 1.283, 0.394, 1.348, 1.445, 0.714, 1.180, 0.603,
                  0.603, 1.545, 0.714, 0.695, 1.236, 1.115, 1.184, 0.808, 1.089, 0.606],
        'default_window': 6,
        'reference': 'Emini EA, Hughes JV, Perlow DS, Boger J. Induction of hepatitis A virus-neutralizing antibody by a virus-specific synthetic peptide. J Virol. 1985;55(3):836-839.',
        'pubmed': '2410642',
        'n_scale': 1.0
    },
    'Parker': {
        'id': 'Parker',
        'title': 'Parker Hydrophilicity',
        'residue': AMINO_ACIDS,
        'score': [2.1, 4.2, 7.0, 10.0, 1.4, 6.0, 7.8, 5.7, 2.1, -8.0,
                  -9.2, 5.7, -4.2, -9.2, 2.1, 6.5, 5.2, -10.0, -1.9, -3.7],
        'default_window': 7,
        'reference': 'Parker JM, Guo D, Hodges RS. New hydrophilicity scale derived from high-performance liquid chromatography peptide retention data. Biochemistry. 1986;25(19):5425-5432.',
        'pubmed': '3539163',
        'n_scale': 1.0
    },
    'Kolaskar-Tongaonkar': {
        'id': 'Kolaskar-Tongaonkar',
        'title': 'Kolaskar and Tongaonkar Antigenicity',
        'residue': AMINO_ACIDS,
        'score': [1.064, 0.873, 0.776, 0.866, 1.412, 1.015, 0.851, 0.874, 1.105, 1.152,
                  1.250, 0.930, 0.826, 1.091, 1.064, 1.012, 0.909, 0.893, 1.161, 1.383],
        'default_window': 7,
        'reference': 'Kolaskar AS, Tongaonkar PC. A semi-empirical method for prediction of antigenic determinants on protein antigens. FEBS Lett. 1990;276(1-2):172-174.',
        'pubmed': '1702393',
        'n_scale': 1.0
    },
    'Chou-Fasman': {
        'id': 'Chou-Fasman',
        'title': 'Chou and Fasman Beta-Turn Prediction',
        'residue': AMINO_ACIDS,
        'score': [0.66, 0.95, 1.56, 1.46, 1.19, 0.98, 0.74, 1.56, 0.95, 0.47,
                  0.59, 1.01, 0.60, 0.60, 1.52, 1.43, 0.96, 0.96, 1.14, 0.50],
        'default_window': 7,
        'reference': 'Chou PY, Fasman GD. Prediction of the secondary structure of proteins from their amino acid sequence. Adv Enzymol Relat Areas Mol Biol. 1978;47:45-148.',
        'pubmed': '364941',
        'n_scale': 1.0
    },
    'Karplus-Schulz': {
        'id': 'Karplus-Schulz',
        'title': 'Karplus and Schulz Flexibility',
        'residue': AMINO_ACIDS,
        'score': [1.046, 1.008, 1.036, 1.068, 0.906, 1.037, 1.094, 1.061, 0.950, 0.927,
                  0.935, 1.102, 0.893, 0.915, 1.049, 1.046, 0.997, 0.904, 0.929, 0.912],
        'default_window': 7,
        'reference': 'Karplus PA, Schulz GE. Prediction of chain flexibility in proteins. Naturwissenschaften. 1985;72:212-213.',
        'pubmed': None,
        'n_scale': 1.0
    },
    'Bepipred': {
        'id': 'Bepipred',
        'title': 'BepiPred Linear Epitope Prediction',
        'residue': AMINO_ACIDS,
        'score': [0.0] * 20,
        'default_window': 0,
        'reference': 'Larsen JE, Lund O, Nielsen M. Improved method for predicting linear B-cell epitopes. Immunome Res. 2006;2:2.',
        'pubmed': '16635264',
        'n_scale': 1.0
    },
    'Bepipred2': {
        'id': 'Bepipred2',
        'title': 'BepiPred 2.0 Linear Epitope Prediction',
        'residue': AMINO_ACIDS,
        'score': [0.0] * 20,
        'default_window': 0,
        'reference': 'Jespersen MC, Peters B, Nielsen M, Marcatili P. BepiPred-2.0: improving sequence-based B-cell epitope prediction using conformational epitopes. Nucleic Acids Res. 2017;45(W1):W24-W29.',
        'pubmed': '28472356',
        'n_scale': 1.0
    }
}

if __name__ == '__main__':
    output_path = os.path.join(os.path.dirname(__file__), 'method_scales.pickle')
    with open(output_path, 'wb') as f:
        pickle.dump(SCALES, f)
    print(f"Created {output_path}")
