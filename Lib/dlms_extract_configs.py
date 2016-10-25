#!/usr/bin/env python
"""Usage:
    dlms_extract_configs.py <nf> [--levcfg <levcfg> --nafion --histpath <hp>]

Extact config frames from DL_MESO HISTORY files, 
convert them to xyz files and place into Dump/ directory.
BEFORE: place DL_MESO script history_config.exe INTO Lib/Extern directory.

Arguments:
    <nf>             Number of config frames requested [default: 10]

Options:
    --levcfg <lc>    Configuration level [default: 0]
    --nafion         Bead ordering "ABCWEP"
    --histpath <hp>  Path of the DL_MESO 'history_config.exe' file

21/06/16
"""
import os, sys
import subprocess
from subprocess import Popen, PIPE
from docopt import docopt


def get_timesteps():
    """Extract number of frames stored in HISTORY binaries
    from CONTROL file."""
    if not os.path.isfile("CONTROL"):
        sys.exit("No CONTROL file present.")
    with open("CONTROL") as f:
        for line in f:
            line = line.rstrip()
            if "traj" in line:
                Nstart, Nevery = list(map(int, line.split()[1:3]))
            if "steps" in line and not "equilibration" in line:
                Nsteps = int(line.split()[1])
    return (Nsteps - Nstart) // Nevery + 1


args = docopt(__doc__)
Nf = int(args["<nf>"])
levcfg = int(args["--levcfg"])

# get number of cores
if not (os.path.isfile("HISTORY") or os.path.isfile("HISTORY000001")):
    sys.exit("No HISTORY* files in the current directory.")
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

if not args["--histpath"]:
    hist_path = "~/Res/PhDcode/Lib/Extern/"
else:
    hist_path = args["--histpath"]
hist_script = os.path.expanduser("~/Res/PhDcode/Lib/Extern/history_config.exe")
transf_script = os.path.expanduser("~/Res/PhDcode/Lib/config2xyz.py")
nafion_flag = "--nafion" if args["--nafion"] else ""
if not os.path.isfile(hist_script):
    sys.exit("Cannot find history_config.exe.")
if not os.path.isfile(transf_script):
    sys.exit("Cannot find config2xyz.py.")

for i in range(Nhf-Nf+1, Nhf+1):
    cmd = "%s %i %i %i" % (hist_script, Nc, levcfg, i)
    print("Running: %s" % cmd)
    subprocess.call(cmd, shell=True)
    subprocess.call("%s CONFIG.out --shift %s" % \
                   (transf_script, nafion_flag), shell=True)
    subprocess.call("mv CONFIG.xyz Dump/dump_%i.xyz" % i, shell=True)


