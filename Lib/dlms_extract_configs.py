#!/usr/bin/env python
"""Usage:
    dlms_extract_configs.py <nf> [--levcfg <levcfg> --nafion --histpath <hp>]

Extact config frames from DL_MESO HISTORY files, 
convert them to xyz files and place into Dump/ directory.
Need to have CONTROL file in the same directory.

Arguments:
    <nf>             Number of config frames requested [default: 10]

Options:
    --levcfg <lc>    Configuration level [default: 0]
    --nafion         Bead ordering "ABCWEP"
    --histpath <hp>  Path of the DL_MESO 'history_config.exe' file

pv278@cam.ac.uk, 21/06/16
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


def get_cores():
    if not (os.path.isfile("HISTORY") or os.path.isfile("HISTORY000000")):
        sys.exit("No HISTORY files in the current directory.")
    pipe = Popen("ls -l HISTORY* | wc -l", stdout=PIPE, shell=True)
    return int(pipe.communicate()[0])


args = docopt(__doc__)
Nf = int(args["<nf>"])
levcfg = int(args["--levcfg"])

# number of frames
Nhf = get_timesteps()
Nc = get_cores()
print("Cores: %i | Total number of frames: %i | Requested frames: %i" % \
      (Nc, Nhf, Nf))
            
if not os.path.exists("Dump"):
    os.makedirs("Dump")

hist_path = "~/sw/bin" if not args["--histpath"] else args["--histpath"]

hist_script = os.path.expanduser(hist_path + "/history_config.exe")
transf_script = os.path.dirname(__file__) + "/config2xyz.py"
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
    subprocess.call("mv CONFIG.xyz Dump/dump_%03i.xyz" % i, shell=True)


