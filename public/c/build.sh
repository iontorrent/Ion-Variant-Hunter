#!/bin/bash

## Find directory of build.sh
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR/src/org/iontorrent/bayesian-scorer
make

mv bayesian-vh-rescorer $DIR
