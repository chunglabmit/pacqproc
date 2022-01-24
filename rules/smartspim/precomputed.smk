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