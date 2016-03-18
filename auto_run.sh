#!/bin/bash

files_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd ..

init_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for temp in "0300" "0600" "0800" "1000" "1200"
do
    dir1=$temp"k"
    mkdir $dir1
    cd $dir1
    for dir2 in 01 02 05 10 15 20
    do
        cp -r $files_dir/template $dir2
        cd $dir2
        sol_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
        for dir3 in 01 05 10 20
        do
            cd $dir3
            
            compa="1.0"
            case "$dir2" in
                01)
                    compa="0.99"
                    ;;
                02)
                    compa="0.98"
                    ;;
                05)
                    compa="0.95"
                    ;;
                10)
                    compa="0.90"
                    ;;
                15)
                    compa="0.85"
                    ;;
                20)
                    compa="0.80"
                    ;;
                 *)
                    echo "solute concentration not listed"
                    exit 1
                    ;;
            esac

            $files_dir/auto_chpar.py kmc_par.h par_temp $temp par_compA $compa
            diff kmc_par.h TEMP_kmc_par.h
            
            echo "in temp:$dir1, sol:$dir2, vcc:$dir3, is it diff OK?(y/n)"
            read yorn
            if [ $yorn != "y" ]; then
                echo "terminating the program ..."
                exit 1
            fi

            cd $sol_dir
        done
    done
    cd $init_dir
done
