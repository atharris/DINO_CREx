# README #

Welcome to the official repository of the Deep-space Interplanetary Navigation Operations Colorado Research Explorer (DINO C-REx)! DINO was developed by a group of graduate researchers through CU Boulder's Graduate Projects course. 

Its purpose is to provide a general-use astrodynamics simulation tool for analyzing optical navigation. To this end, DINO consists of four primary components:
* A Basilisk-based physics engine
* Camera and lighting simulation 
* OpenCV-based image processing and object identification
* Custom batch orbit determination estimators

### How do I get set up? ###

First, ensure that you have the following Python dependencies (available from pip using the command "pip install --user XXXX"):

* pandas
* cv2
* numpy
* scipy
* datetime
* matplotlib

Once these are installed, clone the repo into a folder adjacent to your Basilisk installation.

DINO C-REx consists only of Python modules which are pre-configured in the filesystem.

To start using DINO C-REx, take a peek at either the /dinoScenarios/ directory or the various other DINO subfolders, which contain unit tests and examples to get you started. Within /dinoScenarios/, DINO_main is configured to run simple, single-sim scenarios to help you familiarize with the system. DINO_modeScenarios runs a full, integrated simulation which shifts from propagation to observation and orbit determination as you configure it, and is suitable for general use.

### Contribution guidelines ###

In general, code contributed to DINO C-REx should conform to the PEP-8 standards. The only exception to this is the use of camelCase naming for all variables to maintain consistency with Basilisk. 

### Who do I talk to? ###

Until January 15th 2017, the current point of contact for the maintenence of this repository is Andrew Harris (Andrew.T.Harris@colorado.edu). At that time, the new repo point of contact will be established with the new DINO C-REx team.
