#!/usr/bin/env sh

current_dir=$pwd
cd ../
class_dir=$pwd
cd -

files=$(/bin/ls -1 | grep .C)

for f in $files; do
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"
linenumber=$(eval "sed -n '/MYCODES.h/=' $f") # search line with pathern
if [ -z "$linenuber" ]; then
    continue
sed -i "$linenumber d" $f
sed -i "$linenumber i #include \"$class_dir/MYCODES.h\"" $f
done
