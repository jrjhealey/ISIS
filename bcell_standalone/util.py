"""


"""

import re


def show_method_list():
    print("""
* All available method names:
-----------------------------
1) Chou-Fasman
2) Emini
3) Karplus-Schulz
4) Kolaskar-Tongaonkar
5) Parker
6) Bepipred-1.0
7) BepiPred-2.0
""")


def show_detail_help():
    print("""
detail help here...
""")


def get_sequence(id):
    from urllib.request import urlopen
    if re.search(r'^[A-Z][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9]', id):
        url = "http://www.uniprot.org/uniprot/%s.fasta" % id
        fasta = urlopen(url).read()
        input_sequences = fasta.split('\n')
    else:
        print("* Please check the input ID.")
        exit(0)
    return ''.join(input_sequences[1:])


def print_chart_table(results=None):
    for(sequence, sequence_name), scores in list(results.items()):
        print(("input: {}".format(sequence_name)))
        if scores["epitopes"] and len(scores["epitopes"]) > 0:
            print("Predicted peptides\nNo\tStart\tEnd\tPeptipe\tLength")
            for epitope in scores.get("epitopes"):
                print(("\t".join(map(str, epitope))))
        for row in scores.get("prediction_result"):
            print(("\t".join(row)))


def id_generator(size=6):
    from random import choice
    from string import digits, ascii_lowercase
    chars = digits + ascii_lowercase
    return "".join([choice(chars) for i in range(size)])


def generate_plot(method_name=None, results=None, plot_path=None):
    import numpy as np
    import matplotlib.cbook
    matplotlib.use('Agg')
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure

    filename_suffixes = []

    for i, scores in enumerate(results.values()):
        # get input sequence name which is the second element in the params
        sequence_name = list(results.keys())[i][1]

        fig = Figure(figsize=(6.5, 4.0), facecolor="#F5F5F5")
        ax = fig.add_subplot(1, 1, 1)

        x = []
        y = []
        threshold = scores.get("threshold")
        for row in scores.get("prediction_result")[1:]:
            x.append(int(row[0]))
            if method_name == "Bepipred":
                y.append(float(row[-2]))
            else:
                y.append(float(row[-1]))

        x = np.array(x)
        y = np.array(y)

        ax.plot(x, y)
        ax.set_xlim(xmin=0)

        ax.set_title("{}".format(sequence_name), fontsize=11)
        ax.set_xlabel("Position", fontsize=9)
        ax.set_ylabel("Score", fontsize=9)

        ax.xaxis.set_tick_params(labelsize=9)
        ax.yaxis.set_tick_params(labelsize=9)
        ax.grid(linestyle="dotted")

        ax.fill_between(x, y, np.float(threshold), where=y > np.float(threshold), color='#FFFF00', interpolate=True)
        ax.fill_between(x, y, np.float(threshold), where=y < np.float(threshold), color='#00CC00', interpolate=True)

        ax.axhline(y=threshold, color="r", linewidth=1, label="Threshold")
        ax.legend(loc='upper right', prop={'size': 9})

        random_id = id_generator()
        canvas = FigureCanvas(fig)
        canvas.print_figure("{}/plot_{}.png".format(plot_path, random_id))

        filename_suffixes.append(random_id)

    print(("* A plot has been generated in '{}' directory with 'plot_' prefixed.".format(plot_path)))
    return filename_suffixes
