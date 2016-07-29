#!/usr/bin/env python
"""Usage:
    extract_configs.py [--nfiles <nf> --levcfg <levcfg>]

Options:
    --nfiles <nf>       Number of config files requested [default: 10]
    --levcfg <levcfg>   Configuration level [default: 0]

21/06/16
"""
import os, sys
import subprocess
from subprocess import Popen, PIPE
from docopt import docopt


def get_timesteps():
    """Number of frames stored in HISTORY"""
    if not os.path.isfile("CONTROL"):
        print("No CONTROL file present.")
        sys.exit()
    with open("CONTROL") as f:
        for line in f:
            line = line.rstrip()
            if "traj" in line:
                Nevery = int(line.split()[2])
            if "steps" in line and not "equilibration" in line:
                Nsteps = int(line.split()[1])
    return Nsteps//Nevery + 1


args = docopt(__doc__)
Nf = int(args["--nfiles"])
levcfg = int(args["--levcfg"])

# number of cores
pipe = Popen("ls -l HISTORY* | wc -l", stdout=PIPE, shell=True)
Nc = int(pipe.communicate()[0])

# number of frames
Nhf = get_timesteps()
print("Cores: %i | Total number of frames: %i | Requested frames: %i" % \
      (Nc, Nhf, Nf))
            
if not os.path.exists("Dump"):
    os.makedirs("Dump")

# call history_config.exe Nf times
# reproduce this:
# for i in {292..301}
#     do history_config.exe 20 0 $i; 
#     config2xyz.py CONFIG.out --shift --nafion
#     mv CONFIG.xyz Dump/dump_${i}.xyz
# done

for i in range(Nhf-Nf+1, Nhf+1):
    subprocess.call("/home/pv278/PTFEsim/Scripts/DLMSbin/history_config.exe %i %i %i" %\
                   (Nc, levcfg, i), shell=True)
    subprocess.call("config2xyz.py CONFIG.out --shift --nafion", shell=True)
    subprocess.call("mv CONFIG.xyz Dump/dump_%i.xyz" % i, shell=True)


