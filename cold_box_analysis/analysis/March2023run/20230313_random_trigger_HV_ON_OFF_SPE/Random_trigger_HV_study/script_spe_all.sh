#!/bin/bash

set -e # if something fails, crash everthing and dont keep going

read -p "Are you sure (y/n)? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    return
fi

for subfolder in $(eval "/bin/ls -d  */"); do
    cd $subfolder
    for subruns in $(eval "/bin/ls -d  */ | grep run"); do
        # echo $x$subruns
        cd $subruns
        echo $(pwd)
        eval "root -l -b -q ../../HV_ON_start/run0_argon2x2/giveMeSphe_argon2x2.C+"
        cd ../
    done
    cd ../
done

set +e # now back to normal
