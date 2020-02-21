import os
import sys
import logging
import argparse

from time import sleep

import chimera
from chimera import openModels
from chimera import Molecule
from chimera import runCommand as rc


logging.basicConfig(format="[%(asctime)s] %(levelname)-8s->  %(message)s",
                    level=logging.NOTSET, datefmt='%d/%m/%Y %H:%M:%S %p')
logger = logging.getLogger(__name__)
logging.getLogger().setLevel(logging.INFO)

def sanitise_attname(attribute):
    if "-" in attribute:
        attribute = attribute.replace("-", "")
    if any(char.isdigit() for char in attribute):
        raise ValueError("Digits are not allowed in attribute definitions.")
    # Ensure the first character is lowercase
    if attribute[0].isupper():
        attribute = attribute[:1].lower() + attribute[1:]
    attribute += "Score"
    return attribute

def get_args():
    """Parse command line arguments"""
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-v",
            "--verbose",
            default="1",
            choices=["0", "1", "2"],
            help="Change the verbosity of output (info/debug/errors)."
        )
        parser.add_argument(
            "-i",
            "--inpath",
            action="store",
            help="The file or directory of files to analyse"
        )
        parser.add_argument(
            "-r",
            "--render",
            action="store_true",
            help="Trigger automatic structure colour rendering/cycling."
        )
        parser.add_argument(
            "-p",
            "--profile",
            action="store_true",
            help="Enable script profiling."
        )
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()

def main(args):
    levels_dict = {"0": logging.WARNING, "1": logging.INFO, "2": logging.DEBUG}
    logging.getLogger().setLevel(levels_dict[args.verbose])
    logger.info("Launching {}...".format(__file__))

    logger.info("Reading file(s)...")
    all_data = []
    try:
        with open(args.inpath) as fh:
            for line in fh:
            # do some basic checking that the file is of the correct format
                if not isinstance(eval(line), dict):
                    continue
                all_data.append(eval(line))
    except IOError as err:
        logger.info("No file detected, assuming directory instead.")
        try:
            for filename in os.listdir(args.inpath):
                with open(os.path.join(args.inpath, filename), 'r') as fh:
                    for line in fh:
                    # do some basic checking that the file is of the correct format
                        if not isinstance(eval(line), dict):
                            continue
                        all_data.append(eval(line))
        except OSError as err:
            logger.error("{}\nNo file or directory of files provided/detected.".format(err))
    logger.info("Got epitopes:")
    if logger.getEffectiveLevel() < 30:
        for i in all_data:
            print(logger.info(i))

    models = openModels.list(modelTypes=[Molecule])

    for m in models:
        if len(set(r.id.chainId for r in m.residues)) > 1:
            logger.info("Multichain models detected. Splitting models...")
            rc("split")
            models = openModels.list(modelTypes=[Molecule])


    attributes = []
    for model in models:
        logger.info("Operating on model: {}.{} {}".format(model.id, model.subid, model.name))
        for epitopes in all_data:
            for method in epitopes:
                if (not epitopes[method]) or (epitopes[method] is None):
                    logger.info("No epitopes predicted for {}".format(method))
                    continue
                else:
                    for entry in epitopes[method]:
                        try:
                            index = str(model.sequences()[0]).index(entry[2]) + 1
                            logger.info("Peptide {} found at index {} in model {}.{} {}".format(entry[2], index, model.id, model.subid, model.name))
                        except ValueError:
                            continue

                        attribute = sanitise_attname(method)
                        attributes.append(attribute)
                        for i in range(index, index+len(entry[2])):
                            res = model.residues[i-1]
                            logger.debug("Setting residue {} to {} value {}".format(res, attribute, entry[3]))
                            setattr(res, attribute, float(entry[3]))

# Add some logic to print out the number of epitopes predicted?

# Might be useful to calculate the ranges for the attributes so that useful
# thresholds can be passed to rangecol ahead of time
    if args.render:
        for attr in reversed(set(attributes)):
            rc("rangecol {} min yellow max red novalue white".format(attr))
            sleep(0.5)

    logger.info("Completed successfully.")

if __name__ == "__main__":
    # import cProfile
    # cProfile.runctx('main()', None, locals())
    args = get_args()
    if args.profile:
        from pyinstrument import Profiler
        profiler = Profiler()
        profiler.start()
        main(args)
        profiler.stop()
        print(profiler.output_text())
    else:
        main(args)
