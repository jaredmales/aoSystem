# aoSystem

A program to analyze an astronomical adaptive optics system, both using spatial PSDs and also applying control laws to temporal PSDs.  Primarily focused on post-coronagraph contrast, it can also be used to develop error budgets and Strehl ratio predictions.

See [Males & Guyon (2018)] (https://ui.adsabs.harvard.edu/abs/2018JATIS...4a9001M/abstract) for a complete explanation of the techniques implemented.  Please cite that paper if you use this code in your work!

This program compiles to a stand alone command line program, which can be configured through command line arguments or with configuration files.  You do not need to actually write c++ code to use aoSystem.  One should be able to interact with it from other languages/environments, e.g. python.

## Installation:

This depends on "mxlib", my library of c++ code.  I am constantly changing/updating/extending mxlib for other projects, so I recommend checking out the "aoSystem" branch like so:

```
$ git clone https://github.com/jaredmales/mxlib.git -b aoSystem --single-branch ./
```

Run that command in the directory you want as the top-level directory of the mxlib source repository on your local machine.  Then follow the remaining [installation instructions](https://jaredmales.github.io/mxlib-doc/group__installation.html).

Note: I have installed this on CentOS 6 and 7, and Ubuntu 16 through 20.04.  Joseph Long has successfully compiled it on macOS Darwin (and contributed a bunch of changes to the makefiles so it works) **Update Summer 2021**: mxlib is building on macs again.


## Build:

Once you have setup the mxlib build system (which just means setting up environment variables), you can build aoSystem with:

```
$ make
```

## Run:

See the [User's Guide](https://github.com/jaredmales/aoSystem/blob/master/doc/UserGuide.md).  

Using the config files in the `examples` directory, you can run the program with

```
$ ./aoSystem -c examples/guyon05.conf  
```

To get online help and see all the options:
```
$ ./aosystem -h
```


