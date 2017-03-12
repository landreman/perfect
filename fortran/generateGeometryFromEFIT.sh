#!/bin/bash

# Wrapper script for Matlab geometry generation script so that an executable is in the same file as the perfect executable and hence is likely to be added to $PATH

# Get name of directory where script is located so that we can make relative paths
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE"  ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE"  )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /*  ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE"  )" && pwd )"

echo $PWD
matlab -nodesktop -nosplash -r "usePath='$PWD';run('$DIR/tools/matlab/generateGeometryFromEFIT_script')"

exit 0
