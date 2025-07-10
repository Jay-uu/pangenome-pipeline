#!/usr/bin/env nextflow
def concat_subsamp_samples(sample_ch, proj_name) {
    sample_ch.collectFile(name: "${proj_name}.subsampled.samples", newLine: true, storeDir: "${proj_name}/subsamples")
}
