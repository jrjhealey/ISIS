import os
import sys
import logging
import chimera

from time import sleep

from chimera import openModels
from chimera import Molecule
from chimera import runCommand as rc


logging.basicConfig(format="[%(asctime)s] %(levelname)-8s->  %(message)s",
                    level=logging.NOTSET, datefmt='%d/%m/%Y %H:%M:%S %p')
logger = logging.getLogger(__name__)
logging.getLogger().setLevel(logging.INFO)
logger.info("Launching {}...".format(__file__))

logger.info("Splitting models...")
rc("split")

logger.info("Reading file(s)...")
all_data = []
try:
    with open(sys.argv[1]) as fh:
        for line in fh:
        # do some basic checking that the file is of the correct format
            if not isinstance(eval(line), dict):
                continue
            all_data.append(eval(line))
except IOError as err:
    logger.info("No file detected, assuming directory instead.")
    try:
        for filename in os.listdir(sys.argv[1]):
            with open(os.path.join(sys.argv[1], filename), 'r') as fh:
                for line in fh:
                # do some basic checking that the file is of the correct format
                    if not isinstance(eval(line), dict):
                        continue
                    all_data.append(eval(line))
    except OSError as err:
        logger.error("{}\nNo file or directory of files provided/detected.".format(err))
logger.info("Got epitopes:")
for i in all_data:
    print(i)


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

attributes = []
for model in openModels.list(modelTypes=[Molecule]):
    logger.info("Operating on model: {}.{} {}".format(model.id, model.subid, model.name))
    for epitopes in all_data:
        for method in epitopes:
            if epitopes[method] is None:
                logger.info("No epitopes predicted for {}".format(method))
                continue
            else:
                for entry in epitopes[method]:
                    try:
                        index = str(model.sequences()[0]).index(entry[2]) + 1
                        logger.info("Peptide {} found at index {}".format(entry[2], index))
                    except ValueError:
                        continue

                    attribute = sanitise_attname(method)
                    attributes.append(attribute)
                    for i in range(index, index+len(entry[2])):
                        res = model.residues[i-1]
                        logger.info("Setting residue {} to {} value {}".format(res, attribute, entry[3]))
                        setattr(res, attribute, float(entry[3]))

for attr in set(attributes):
    rc("rangecol {} min yellow max red novalue white".format(attr))
    sleep(0.5)
