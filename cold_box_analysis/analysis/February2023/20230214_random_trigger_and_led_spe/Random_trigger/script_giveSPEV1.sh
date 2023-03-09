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
    if [[ $subruns == *"run99"*  ||  $subruns == *"run0"*   ||  $subruns == *"run1"* ]]; then
        continue
    fi
    cd $subruns
    echo $(pwd)
    eval "root -l -b -q ../run2_switching_on_argon2x2_argon4_DCem1dot2VD_47V_DCem1dot0VD_36V/giveMeSpheV1.C+"
    cd ../
done

set +e # now back to normal
