#!/bin/bash

## Note, sbcl has to be installed into the system.
## 
## sudo apt-get install sbcl slime emacs on Ubuntu
## or sbcl.org for other platforms
##
## Note, if you use the Ubuntu specific package, the compiled binary 
## may not work in other LINUX like CentOS

## First argument is optional and can specify the name of the binary to build

echo
echo Building ion-variant-hunter-core from src. . . 
echo
echo Note, this build script requires that sbcl is installed.
echo If on ubuntu, can do \'sudo apt-get install sbcl slime\', note however
echo if done on 10.04, the resulting binary won\'t run on CentOS 5.  If you
echo would like to have it run on both platforms, download the tarball of
echo sbcl from www.sbcl.org and install it into your system.
echo 

## Find directory of build.sh
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR/src

sbcl --eval '(progn (load "system.lisp") (build "'$1'"))'

mv ion-variant-hunter-core ..
