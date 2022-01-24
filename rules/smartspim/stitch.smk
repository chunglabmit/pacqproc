rule stitch_compute_master:
    input: f"{config['master_channel']}{destripe_flag}/.done"
    output:
        import_xml=f"{config['master_channel']}{destripe_flag}/import.xml",
        displacement_xml=f"{config['master_channel']}{destripe_flag}/displacement.xml",
    threads: nthreads
    params: subvoldim=200
    run: 
        master_path = Path(config['workdir'])/f"{config['master_channel']}{destripe_flag}"
        shell("{terastitcher} -1 --volin={master_path} \
                                 --ref1=H --ref2=V --ref3=D \
                                 --vxl1={vx} --vxl2={vy} --vxl3={vz} \
                                 --projout={output.import_xml} --sparse_data")

        mpstitcher = config['mpstitcher']
        os.environ['N_PROCESSORS'] = str(threads)
        print(threads, 'threads')
        shell("python {mpstitcher} -2 --projin={output.import_xml} --projout={output.displacement_xml} --subvoldim={params.subvoldim}")
        shell("rm -f *.out") # clean up intermediate files


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