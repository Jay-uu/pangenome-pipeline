#!/usr/bin/env nextflow
def summarize_bintables(bintable_ch, proj_name) {
    def bintable_header = Channel.value("Bin ID	Completeness	Contamination	Tax GTDB-Tk")
    def bintable_rows = bintable_ch
        .collectFile(keepHeader: true, skip: 2)
        .splitCsv(header: true, skip: 1, sep: "\t")
        .map { row -> "${row.'Bin ID'}	${row.Completeness}	${row.Contamination}	${row.'Tax GTDB-Tk'}" }
    bintable_header.concat(bintable_rows).collectFile(name: "summarized_bins.bintable", newLine: true, sort: false, storeDir: "${proj_name}/bins")
}
