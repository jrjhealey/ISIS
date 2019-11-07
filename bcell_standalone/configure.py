#! /usr/bin/env python
"""
Created on: 05/22/2017

@author: Dorjee Gyaltsen
@brief: This script will download and install bcell-executable package.
"""


import os
import sys
import pip
import subprocess
import pkg_resources
from distutils.spawn import find_executable

NETSURFP_BIN  = os.environ.get('NETSURFP_BIN',  'netsurfp')

class Configure(object):

    def __init__(self, version="2.0"):
        self.version = version
        self.project_dir = os.path.abspath(".")
        self.deps = os.path.join(self.project_dir, "deps")

    def install_requirements(self):
        requirements = os.path.join(self.project_dir, "requirements.txt")
        try:
            print("Installing requirements file (this might take some time)...")
            subprocess.check_call(["pip", "install", "-r", requirements])
        except subprocess.CalledProcessError as e:
            print(e.output)

    def check_prerequisites(self):
        prerequisites = ["tcsh", "gawk"]
        for prereq in prerequisites:
            
            executable = find_executable(prereq)
            if not executable:                
                print("* prerequisite: '{}' is not found in your path.  Please install and then rerun the configure script".format(prereq))
                sys.exit(1)                
            else:
                print ('prerequisite: "%s" was found.' % prereq)
        # if we get here, all prerequisites have been found
        print("All prerequisites found!")

    def check_netsurfp(self):
        executable = find_executable(NETSURFP_BIN)
        if not executable:
            print('* optional dependency: "netsurfp" was not found in your path. "netsurfp" is a dependency for method "Bepipred-2.0". Please install it as the direction in the "README" file before you run prediction with method "Bepipred-2.0".\n')
        else:
            print ('optional dependency: "%s" was found.\n' % NETSURFP_BIN)

    def check_python_modules(self):
        missing_packages_list = []
        installed_packages = pkg_resources.working_set
        installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
        required_packages_list = ['numpy', 'scipy', 'scikit-learn==0.18', 'matplotlib==2.0.0']
        for package in required_packages_list:
            if not any([p.startswith(package) for p in installed_packages_list]):
                missing_packages_list.append(package)
        return missing_packages_list

if __name__ == "__main__":
    config = Configure()
    config.check_prerequisites()
    config.check_netsurfp()
    # un-comment the line below, if you like configure script to take care of the required package install
    config.install_requirements()
    missing_packages_list = config.check_python_modules()
    if missing_packages_list:
        print("* You must have %s package(s) installed.\nRun this command to install them:\n$ pip install %s \n" % (', '.join(missing_packages_list), ' '.join(missing_packages_list)))
    else:
        print("That's it. You're all set!")

