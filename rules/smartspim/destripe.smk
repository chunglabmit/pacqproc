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