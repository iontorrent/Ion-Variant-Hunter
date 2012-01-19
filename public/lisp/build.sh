#!/bin/bash

## Note, sbcl is a requirement to build this system.
## 
## sudo apt-get install sbcl slime emacs on Ubuntu
## or sbcl.org for other platforms

## First argument is optional and can specify the name of the binary to build

echo
echo Building ion-variant-hunter-core from src. . . 
echo
echo This build script requires that sbcl is installed.
echo For greatest stability, download source from www.sbcl.org, and make modifications specified in
echo the sbcl-mod directory.
echo
echo Note, building under Ubuntu 10.04 will result in a binary won\'t run on CentOS 5.
echo

## Find directory of build.sh
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR/src/org/iontorrent/variant-hunter/

sbcl --eval '(progn (load "system.lisp") (build "'$1'"))'

mv $DIR/src/org/iontorrent/variant-hunter/ion-variant-hunter-core $DIR
cp -p $DIR/scripts/samRegionOverlap.py $DIR
