#!/bin/bash

## Find directory of build.sh
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIN=bayesian-vh-rescorer

rm $DIR/$BIN
