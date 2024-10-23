#! /bin/bash

N=1000
Q=0.0
SEED=0


source venv/bin/activate

./PlummerPlus/PlummerPlus.py -n ${N} -q ${Q} -rs ${SEED}

julia ConvertToHenon.jl --N ${N} --q ${Q} --seed ${SEED}

rm output.txt

deactivate
