#!/bin/bash

set -e # if something fails, crash everthing and dont keep going

read -p "Are you sure (y/n)? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    return
fi

for subruns in $(eval "/bin/ls -d  */ | grep run"); do
    # echo $x$subruns
    if [[ $subruns == *"run99"*  ||  $subruns == *"run78"*   ||  $subruns == *"run1"* ]]; then
        continue
    fi
    cd $subruns
    echo $(pwd)
    eval "root -l -b -q ../run74_all_devices/giveMeSpheV1.C+"
    cd ../
done

set +e # now back to normal