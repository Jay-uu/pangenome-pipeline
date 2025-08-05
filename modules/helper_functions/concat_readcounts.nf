#!/usr/bin/env nextflow
def concat_readcounts(readcounts_ch, proj_name) {
    readcounts_ch.collectFile(keepHeader: true, name: "original_readcounts.tsv", newLine: false, sort: false, storeDir: "${proj_name}/subsamples")
}
