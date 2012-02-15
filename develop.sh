#!/bin/bash

SRC_LOC=org/iontorrent/variant-hunter
DEV_DIR=public/lisp/src/$SRC_LOC
TEST_DIR=public/lisp/test/$SRC_LOC

echo 'This script builds the project, and then runs emacs loading'
echo 'the .lisp file for development.'
echo
echo 'Now building the projects for required JAVA components'
echo '******************************************************'
ant
cd dist
echo
echo '******************************************************'
echo "Now loading up emacs on X display $DISPLAY"
emacs ../$DEV_DIR/*.lisp ../$TEST_DIR/*.lisp ../public/lisp/scripts/* &> /dev/null &> /dev/null &
echo
echo 'Emacs hopefully launched now.  If SLIME is also installed,'
echo 'one can start up an interactive LISP shell, but doing,'
echo 'Meta-x slime'
echo
