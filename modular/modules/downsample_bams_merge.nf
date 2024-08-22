/*
This process takes the SqueezeMeta output and contig names to downsamples the bam files to get even coverage between samples, in preparation for Variant Calling.
Input: the results directory from SqueezeMeta.
Output: Since there is a possibility that no bams fit the minimum coverage and breadth criteria, this process might have no output or it will send
        a tuple with a fasta file of all contigs longer than 1000 bases from the input pangenome and a merged bam-file from all bams that passed the breadth
        and coverage criteria.
The downsampling shell code is modified from POGENOM's Input_pogenom pipeline by Anders Andersson and Luis F. Delgado
See here: https://github.com/EnvGen/POGENOM/blob/master/Input_POGENOM/src/cov_bdrth_in_dataset.sh
*/
process downsample_bams_merge {
    publishDir "${params.project}/mOTUs/results/${pang_id}/pangenome", pattern: "cov_breadth.txt"
    label "downsample_bams_merge"
    tag "${pang_id}"
    input:
    tuple(val(pang_id), path(pang_sqm), path(contigs_tsv))
    output:
    tuple(path("${pang_sqm}_long_contigs.fasta"), path("${pang_sqm}_merged.bam"), optional: true, emit: ref_merged)
    path("NOT_PASSED.txt"), optional: true, emit: not_passed_message
    path("cov_breadth.txt"), emit: cov_breadth
    shell:
    '''
    if [ ! -s !{contigs_tsv} ]; then
        #create a contigs tsv will all contigs and bed
        #required fields: name, start, length
        bam=$(ls !{pang_sqm}/data/bam/*.bam | head -n1) #*/ comment is for editor purposes
        echo "Using ${bam} to create contigs.bed, since no selected contigs were provided."
        bedtools bamtobed -i $bam > contigs.bed
    else
        #create bed from contigs_tsv
        echo "Converting provided list of contigs to a bed file."
        awk ' { print $1, 1, $2} ' !{contigs_tsv} > contigs.bed
    fi
        
    #Create tmp bams, filter for contigs over ${cont_len} bases put reads aligning to them in tmp_bams
    echo "Creating tmp bams"
    mkdir tmp_bams
    for bam in !{pang_sqm}/data/bam/*.bam; #*/
    do
	echo "Filtering ${bam} alignments for selected contigs"
        #Filter to select only paired reads (-f 2) and avoids optical duplicates (-F 1024)
        samtools index ${bam}
        bam_ID=$(basename $bam .bam)
        samtools view -Sbh -F 1024 -q 1 -L contigs.bed --threads !{task.cpus} $bam > tmp_bams/${bam_ID}.bam
    done
    
    mkdir -p !{pang_sqm}_mergeable
    echo "Creating depth files and checking breadth and coverage."    
    #create depth files, col 3 is number of reads at one position
    echo -e "Genome\tSample\tMedian_cov\tBreadth\n" > cov_breadth.txt
    for bam in tmp_bams/*.bam;
    do
        #mpileup command doesn't allow multithreading
        #-A for count orphans
        #samtools mpileup -A -d 1000000 -Q 15 -a $bam > tmp.mpileup
        #using samtools depth instead, newer
        samtools depth -aa -J -s -q 1 -b contigs.bed $bam > tmp.depth
    
        # ---- arguments
        #mpileupfile=tmp.mpileup
        outbamfile=$(basename $bam bam)downsampled.bam #name of output
        mag=!{pang_sqm} #pangenome name
        mincov=!{params.min_cov}
        minbreadth=!{params.min_breadth}
        samplename=$(basename ${bam#"${mag}."} .bam) #this might be bugged. When testing manually it's just the samplename but my output file is mag.sample. Not removing string successfully.

        #--- Median coverage
        #col 4 has nr of reads mapped to position, only take positions where reads mapped, sort by numerical value, add to array,
        #take value in middle of array (or mean of two middle values if even nr of values) = median cov
        #Not 100% sure why 0 positions are excluded in Input_pogenom. Their preprint does say that they do it purposefully though.
        #https://www.biorxiv.org/content/10.1101/2020.03.25.999755v1.full
        #but assuming it's because they aren't actually used for variant calling and therefore irrelevant for the coverage and downsampling
        #cov=$(cut -f4 $mpileupfile | grep -vw "0" | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
        #samtools depth has the nr of reads at position in col 3
        cov=$(cut -f3 tmp.depth | grep -vw "0" | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')

        #---breadth
        #non_zero=$(cut -f4 $mpileupfile | grep -cvw "0") #mpileup way
        non_zero=$(cut -f3 tmp.depth | grep -cvw "0") #depth way
        positions=$(wc -l < tmp.depth)
        breadth=$(echo $non_zero*100/$positions | bc -l )

        echo "Genome:" $mag "- Sample:" $samplename "Median_coverage of core:" $cov " breadth %:" $breadth
        echo -e "${mag}\t${samplename}\t${cov}\t${breadth}\n" >> cov_breadth.txt

        #---selection of BAM files and downsample
        if (( $(echo "$breadth >= $minbreadth" | bc -l) )) && (( $(echo "$cov >= $mincov" | bc -l) )); then
            echo "Downsampling coverage to $mincov - Genome: $mag - Sample: $samplename "
            if (( $(echo "$mincov == 0" | bc -l) )); then
                $mincov=$cov
            fi
            limite=$(echo "scale=3; $mincov/$cov" | bc )
            samp=$(echo "scale=3; ($limite)+10" | bc)
            samtools view -Sbh --threads !{task.cpus} -s $samp $bam | samtools sort -o !{pang_sqm}_mergeable/$outbamfile --threads !{task.cpus}
        fi
    done
    
    #Merge bam-files that pass the check, if more than one bam in mergeable/ #*/ what to do if no files?
    #if at least one file in mergeable, create new fasta with only long contigs
    echo "Checking mergeable"
    if [ -z "$(ls -A !{pang_sqm}_mergeable)" ]; then
         echo "No sample fit the alignment criteria. Skipping further analysis for !{pang_sqm}"
         cat "WARNING: No sample fit the alignment criteria for !{pang_sqm}. If you want to analyze this sample further try lowering --min_cov and/or --min_breadth." > NOT_PASSED.txt
    else
        echo "Merging downsampled bams. and creating fasta of pangenome with only NBPs over !{params.min_contig_len} bases."
        ls !{pang_sqm}_mergeable/*.bam > bamlist.txt #*/
        samtools merge -o !{pang_sqm}_merged.bam -b bamlist.txt --threads !{task.cpus}
        samtools index !{pang_sqm}_merged.bam --threads !{task.cpus}
        samtools idxstats !{pang_sqm}_merged.bam --threads !{task.cpus}| awk '$2 >= '!{params.min_contig_len}' { print $0 }' > long_contigs.tsv
        awk '{ print $1 }' long_contigs.tsv > contig_names.tsv
        #seqtk doesn't allow multithreading
        seqtk subseq !{pang_sqm}/results/01.*.fasta contig_names.tsv > !{pang_sqm}_long_contigs.fasta
    fi
    #cleanup step
    rm -r tmp*
    '''

}
