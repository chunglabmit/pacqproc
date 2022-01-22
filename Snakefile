#!/usr/bin/env python3
from pathlib import Path
import os

configfile: "config.yaml"
workdir: config['workdir']

# Variable settings
vx, vy, vz = map(float, config['voxel_size'].split(','))
destripe_flag = "_destriped" if "destripe" in config['processes'] else ""
stitch_flag = destripe_flag + "_stitched"
nthreads = min(os.cpu_count(), int(config['threads']))
os.environ['N_PROCESSORS'] = str(nthreads)
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


rule destripe:
    output: touch("{channel}_destriped/.done")
    params:
        sigma1 = 256,
        sigma2 = 256,
        wavelet = "db5",
        crossover = 10
    threads: nthreads
    run:
        shell("pystripe --input {wildcards.channel} \
                --sigma1 {params.sigma1} --sigma2 {params.sigma2} \
                --wavelet {params.wavelet} --crossover {params.crossover} \
                --workers {threads}")

        # Prepare destriped image tiles to be compatible with Terastitcher
        # 1. fills up missing z-slices
        # 2. make the xyz coords directory names positive numbers, if needed
        shell("tsv-fill-blanks --src {wildcards.channel}")
        shell("tsv-renumber-directories --path {wildcards.channel}_destriped")


rule stitch_compute_master:
    input: f"{config['master_channel']}" + destripe_flag + "/.done"
    output:
        import_xml=f"{config['master_channel']}" + destripe_flag + "/import.xml",
        displacement_xml=f"{config['master_channel']}" + destripe_flag + "/displacement.xml",
    threads: nthreads
    run: 
        master_path = Path(config['workdir'])/config['master_channel']
        shell("{terastitcher} -1 --volin={master_path} \
                                 --ref1=H --ref2=V --ref3=D \
                                 --vxl1={vx} --vxl2={vy} --vxl3={vz} \
                                 --projout={output.import_xml} --sparse_data")

        mpstitcher = config['mpstitcher']
        shell("python {mpstitcher} -2 --projin={output.import_xml} --projout={output.displacement_xml}")


rule copy_master_displacement:
    input: f"{config['master_channel']}" + destripe_flag + "/displacement.xml",
    output: expand("{channel}" + destripe_flag + "/displacement.xml", channel=config['non_master_channels']),
    run:
        xmlpath = Path(config['workdir'])/input[0]
        for ch in config['non_master_channels']:
            dest = Path(config['workdir'])/(ch + destripe_flag)/xmlpath.name
            shell("cp {xmlpath} {dest}")
            shell("xmlstarlet ed -L -u '//TeraStitcher/stacks_dir/@value' -v {dest.parent} {dest}")


rule stitch_project:
    input: f"{{channel}}{destripe_flag}/displacement.xml"
    output: f"{{channel}}{destripe_flag}/displproj.xml"
    shell: "{terastitcher} -3 --projin={input} --projout={output}"


rule stitch:
    input: f"{{channel}}{destripe_flag}/displproj.xml"
    output: touch(f"{{channel}}{destripe_flag}_stitched/.done")
    threads: nthreads
    params:
        compression = config['compression']
    run:
        output_pattern = Path(output[0]).parent/'img_{z:04d}.tiff'
        shell("tsv-convert-2d-tif --xml-path {input}\
                --output-pattern \"{output_pattern}\" \
                --cpus {threads} --compression {params.compression} \
                --ignore-z-offsets")
        

rule convert_to_precomputed:
    input: expand("{channel}" + stitch_flag + "/.done", channel=config['channels'])
    output: touch(expand("{label}.precomputed/.done", label=config['labels']))
    params:
        levels = 7,
        format = "blockfs"
    threads: min(nthreads, 32)
    run:
        for channel, label in zip(config['channels'], config['labels']):
            input_pattern = channel + stitch_flag + "/img_*.tif*"
            outpath = label + ".precomputed"
            voxel_size = ",".join(map(str, [vx, vy, vz]))
            shell("precomputed-tif \
                    --source \"{input_pattern}\" --dest {outpath} \
                    --levels {params.levels} --format {params.format} \
                    --voxel-size {voxel_size} --n-cores {threads}")