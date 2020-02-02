#! /usr/bin/env python
"""
Created on: 05/03/2017

@author: Dorjee Gyaltsen
@brief: This script will download and install population-coverage-pickle package.
"""
from __future__ import print_function

import os
import subprocess


class Configure(object):

    def __init__(self, version="2.1.3"):
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

if __name__ == "__main__":
    config = Configure()

    # un-comment the line below, if you like configure script to take care of the required package install
    config.install_requirements()

    print("* You must have 'numpy' and 'matplotlib' packages installed.\nRun this command to install them:\n"
    "$ pip install numpy matplotlib==2.0.0\n")

    print("That's it. You're all set!")
