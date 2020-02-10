#!/usr/bin/env python2.7
"""
A wrapper script to run the various tools within the
in silico immunogenicity IEDB suite
"""
from __future__ import print_function
import argparse
import logging
import sys

from bcell_standalone.BCellRunner import BCellRunner

logging.basicConfig(format="[%(asctime)s] %(levelname)-8s->  %(message)s",
                    level=logging.NOTSET, datefmt='%d/%m/%Y %H:%M:%S %p')
logger = logging.getLogger(__name__)


logging.getLogger().setLevel(logging.INFO)
logger.info("Launching {}...".format(__file__))

# Run the BCellStandalone tools
logger.info("Calling BCellRunner...")

BCR = BCellRunner(sys.argv[1], sys.argv[2])
epitopes = {method:BCR.results[method].values()[0]['epitopes'] for method in BCR.results}
print(epitopes)
