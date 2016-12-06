#!/bin/bash

for i in */; do sed -i s/sourcesVStructure\ =\ 2\ 1/sourcesVStructure\ =\ 1\ 2/ $i/input.namelist; done

for i in */; do sed -i s/gConstraints\ =\ 2\ 1/gConstraints\ =\ 1\ 2/ $i/input.namelist; done

for i in */; do sed -i s/solverTolerance\ =\ 1d-[0-9]*/solverTolerance\ =\ 1d-10/ $i/input.namelist; done

