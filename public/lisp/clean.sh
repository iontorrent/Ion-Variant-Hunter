#!/bin/bash

## Find directory of build.sh
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

BIN=ion-variant-hunter-core

rm $DIR/$BIN $DIR/samRegionOverlap.py
