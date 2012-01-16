#!/bin/bash

## Find directory of build.sh
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR/bayesian-scorer
make

mv bayesian-vh-rescorer ..
