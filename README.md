[![Homepage](https://img.shields.io/badge/Home-plumed.org-green.svg)](http://www.plumed.org)
[![Homepage](https://img.shields.io/badge/Google_group-plumed--users-green.svg)](http://groups.google.com/forum/#!forum/plumed-users)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](http://www.gnu.org/licenses/lgpl-3.0)
[![Github Releases](https://img.shields.io/github/release/plumed/plumed2.svg)](https://github.com/plumed/plumed2/releases)
[![Twitter Follow](https://img.shields.io/twitter/follow/plumed_org.svg?style=social&label=Follow)](https://twitter.com/plumed_org)

Branches and releases
---------------------

Several branches and tags are stored on the git repository.

Branches named `v2.X` correspond to release branches.

Master branch may contain non tested features and is not expected to be used by non-developers.
It typically contains features that will be available on the next release.

Tags named `v2.XbY` correspond to beta releases, use it with care.
Tags named `v2.X.Y` correspond to official releases, use the latest available.

In addition, the repository contains a number of other branches related to specific features.
Please contact the developers that are committing on those branches before basing your work
there, since they might contain temporary work and might be rebased later.
For instance, branch `testdoc` is setup so as to push a test copy of the manual
and is often force pushed.

To report problems found on beta or official releases, use the normal
[plumed-users@googlegroups.com](mailto:plumed-users@googlegroups.com)
mailing list. Please state exactly which version you are using.
To report problems found on `master` branch, use the
[plumed2-git@googlegroups.com](plumed2-git@googlegroups.com) mailing list.
This is also the correct place for discussions about new features etc.
When reporting please provide the git hash (you can obtain it with `git rev-parse HEAD`).

Status
------

Below you find the status on [Travis-CI](http://travis-ci.org/plumed/plumed2) for the release branches.

| Branch   |      Status   | First stable release | Still supported |
|:--------:|:-------------:|:--------:|:------:|
| master   | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=master)](https://travis-ci.org/plumed/plumed2) | v2.5.0 is expected late 2018 | / |
| v2.4     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.4)](https://travis-ci.org/plumed/plumed2)   | v2.4.0 is expected late 2017 | yes |
| v2.3     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.3)](https://travis-ci.org/plumed/plumed2)   | Dec 12, 2016 | yes |
| v2.2     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.2)](https://travis-ci.org/plumed/plumed2)   | Oct 13, 2015 |no |
| v2.1     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.1)](https://travis-ci.org/plumed/plumed2)   | Sep 9, 2014 | no |
| v2.0     | Not available | Sep 27, 2013 | no |

Content
-------

Here's a description of the content of each file and directory in the root PLUMED directory.

    CHANGES          : change log
    COPYING.LESSER   : license
    Makefile         : makefile
    Makefile.conf.in : template configuration makefile
    PEOPLE           : list of authors
    README           : this file
    VERSION          : version file
    astyle           : a local version of astyle, used to format code
    configurations   : template configuration files
    configure        : configuration script
    configure.ac     : configuration script (autoconf)
    developer-doc    : developer documentation
    include          : symbolic link for include files
    macports         : directory where Portfiles are generated
    patches          : patch scripts
    release.sh       : developer utility to publish releases
    regtest          : regression tests, including reference results
    scripts          : shell tools
    src              : source code
    sourceme.sh.in   : template configuration script
    test             : examples
    user-doc         : user documentation
    vim              : directory where vim syntax is generated

Required software
-----------------

Required software:

* GNU make.
* C/c++ compiler (c++11 support is required as of version 2.4).
* A modern version of the `patch` command line tool.
* Support for POSIX library `dirent.h`.
* `xxd` (present in most UNIX distributions).

Suggested software (libraries are checked by `./configure` and enabled if available):

* MPI library to run parallel simulations. It should be the same library used by your MD code.
* Optimized blas and lapack libraries. They are automatically replaced by an internal version if not available.
* [VMD molfile plugins](http://www.ks.uiuc.edu/Research/vmd/plugins) to read arbitrary file formats. They are automatically replaced by an internal version supporting a few formats if not available.
* [Zlib library](http://zlib.net/) to use compressed data files.
* [Xdrfile library](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library) to have read/write access to gromacs
  trajectory files.
* [Doxygen](http:://www.doxygen.org) to build user manual. Doxygen might need the following packages:
  * Latex to build the pdf user manual.
  * [Graphviz](http://www.graphviz.org) to show class hierarchy in
    developer manual.

Quick compilation instructions
------------------------------

Extensive installation instructions are in the [user documentation](http://www.plumed.org/documentation).
Configure for your system

    ./configure --prefix=$HOME/opt
    
    
If necessary, edit `Makefile.conf`. 
Configure your environment

    source ./sourceme.sh
    
Compile plumed

    make
    
The `plumed` executable should be now in your execution path

    plumed help
    
Compile the manuals.

    make doc

User documentation can be found at `user-doc/html/index.html`.
Developer documentation can be found at `developer-doc/html/index.html`.
[Pre-compiled documentation](http://www.plumed.org/documentation) is available online, so this is only required
if you are working with a modified version of the code!

Install PLUMED in `$HOME/opt` (directory should be set during `./configure`):

    umask 022
    make install
    
A sample modulefile with environment variable will be placed in
`$HOME/opt/lib/plumed/src/lib/modulefile`. If you compiled your own documentation, paths to the installed documentation can be found with command `plumed info --user-doc`.
