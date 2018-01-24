# aoSystem

A program to analyze an astronomical adaptive optics system, both using spatial PSDs and also applying control laws to temporal PSDs.  Primarily focused on post-coronagraph contrast, it can also be used to develop error budgets and Strehl ratio predictions.

See https://arxiv.org/abs/1712.07189 for a complete explanation of the techniques implemented.  Please cite that paper (accepted to JATIS) if you use this code in your work!

This program compiles to a stand alone command line program, which can be configured through command line arguments or with configuration files.  You do not need to actually write c++ code to use aoSystem.  One should be able to interact with it from other languages/environments e.g. python.

## Installation:

This depends on "mxlib", my library of c++ code.  I am constantly changing/updating/extending mxlib for other projects, so I recommend checking out the "aoSystem" branch like so:

git clone https://github.com/jaredmales/mxlib.git -b aoSystem --single-branch ./

Run that command in the directory you want as the top-level directory of mxlib.  Then follow the remaining installation instructions at: https://jaredmales.github.io/mxlib/group__installation.html

Note: I have installed this on CentOS6 and 7, and Ubuntu 16.04 and 17.04.  Joesph Long has successfully compiled it on macOS Darwin (and contributed a bunch of changes to the makefiles so it works).



## Build:

Once you have setup the mxlib build system (which just means setting up environment variables), you can build aoSystem with:

make -B -f $MXMAKEFILEINC t=aoSystem

## Run:

Using the example config file, you can run the program with

> ./aoSystem -c example.conf  

See the documentation in example.conf.  Command line options can be interrogated with

> ./aosystem -h



