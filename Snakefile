#!/usr/bin/env python3
from pathlib import Path
import os

configfile: "config.yaml"
workdir: config['workdir']

# Variable settings
vx, vy, vz = map(float, config['voxel_size'].split(','))
destripe_flag = "_destriped" if "destripe" in config['processes'] else ""
stitch_flag = f"{destripe_flag}_stitched"
nthreads = min(os.cpu_count(), int(config['threads']))
terastitcher = config['terastitcher']
mpstitcher = config['mpstitcher']


# Target settings
targets = []
if "destripe" in config['processes']:
    targets.append(expand("{channel}_destriped/.done", channel=config['channels']))
if "stitch" in config['processes']:
    config['master_channel'] = config['channels'][config['master_channel_index']]
    config['non_master_channels'] = set(config['channels']) - set([config['master_channel']])
    targets.append(expand("{channel}" + destripe_flag + "/displacement.xml", channel=config['channels']))
    targets.append(expand("{channel}" + destripe_flag + "/displproj.xml", channel=config['channels']))
    targets.append(expand("{channel}" + stitch_flag + "/.done", channel=config['channels']))
if "precomputed" in config['processes']:
    targets.append(expand("{label}.precomputed/.done", label=config['labels']))


# Other settings
wildcard_constraints: 
    channel=r"^.+?(?=_destriped)"


# Rules
rule all:
    input: targets

if config['workflow'] == "smartspim":
    if "destripe" in config['processes']:
        include: "rules/smartspim/destripe.smk"
    if "stitch" in config['processes']:
        include: "rules/smartspim/stitch.smk"
    if "precomputed" in config['processes']:
        include: "rules/smartspim/precomputed.smk"