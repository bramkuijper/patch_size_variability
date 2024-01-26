#!/usr/bin/env python3

import numpy as np
import datetime

baseline_mortality = np.linspace(0.01,1.0,num=30)

n = [[5,5]]

d = 0.5

cfec = 0.01

bfec = np.linspace(0.001,0.14,num=30)

exe = "patch_size.exe"


date = datetime.datetime.now()
base_name = "sim_patch_size_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

ctr = 0

for mort_i in baseline_mortality:
    for n_i in n:
        for bfec_i in bfec:

            ctr += 1
            base_name_i = base_name + "_" + str(ctr)

            print(f"./{exe} " +\
                    f"{d} " +\
                    f"{n_i[0]} " +\
                    f"{n_i[1]} " +\
                    f"{mort_i} " +\
                    f"{bfec_i} " +\
                    f"{cfec} " +\
                    f"{base_name_i} "
                    )

