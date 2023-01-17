#!/usr/bin/env sh

current_dir=$( pwd )
cd ../
class_dir=$( pwd )
cd $current_dir

files=$(/bin/ls -1 | grep .C)

for f in $files; do
    linenumber=$(eval "sed -n '/MYCODES.h/=' $f") # search line with pathern
    if [ -z "$linenumber" ]; then
        continue
    fi
    sed -i "$linenumber d" $f
    sed -i "$linenumber i #include \"$class_dir/MYCODES.h\"" $f
done


cd ../
COUNTER=0
lines=$(eval "sed -n '/define/=' load_my_class.sh")

for line in $lines; do
    sed -i "$line d" load_my_class.sh
    if [ $COUNTER -eq 0 ]; then
        sed -i "$line i \ \ \ \ eval \'root -e \"#define memorydepth \$npts\" -e \".L $class_dir/MYCODES.h\"\'" load_my_class.sh
    else
        sed -i "$line i \ \ \ \ eval \'root -e \"#define memorydepth \$npts\" -e \".L $class_dir/MYCODES.h\" -e \"SAMPLE s\(\\\\\"s\\\\\",\\\\\"\$file\\\\\"\)\"\'" load_my_class.sh
    fi
    COUNTER=$(( COUNTER + 1 ))

done
