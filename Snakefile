

rule linear_bcell_prediction:
    input:
        "data/sequences/PAK_01787.faa"

    output:
        "output/PAK_01787.lbcp"
    shell:
        "python2.7 bcell_standalone/predict_antibody_epitope.py -f ../data/sequences/PAK_01787.faa -m Chou-Fasman"
