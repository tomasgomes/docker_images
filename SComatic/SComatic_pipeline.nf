/*
 * Defines necessary parameters
 */
params.samplename = false
//params.samplename = "ID"
params.outdir = "./"
params.genome = false
//params.genome = "refdata-gex-GRCh38-2020-A/fasta/genome.fa"
params.bam = false
//params.bam = "outs/possorted_genome_bam.bam"
params.ct = false
//params.ct = "ct_df.tsv"
params.mincells = 4
params.mincelltypes = 2
params.cores = 12



/*
 * Input parameters validation and preparation
 */
if(!params.bam) exit 1, "Please provide an input sorted BAM file"
if(!params.ct) exit 1, "Please provide an input cell type annotation. Should be a tab-separated table with Index and Cell_type columns"
if(!params.samplename) exit 1, "Please provide a name for this sample"
if(!params.genome) exit 1, "Please provide a genome reference in FASTA format"

if(file("${params.outdir}${params.samplename}").exists()) println "Output folder already exists. No files will be overwritten, but execution may fail."



/*
 * Step 0. Getting things going
 */
process starting{

    script:
    """
    echo Starting...
    """
}

/*
 * Step 1. Split the BAM file by cell type
 */
process splitbam {

    input:
    file(params.bam)
    //file transcriptome_file

    output:
    file params.transindex into transcriptome_index

    script:
    """
    python SComatic/scripts/SplitBam/SplitBamCellTypes.py --bam ${params.bam} \\
        --meta ${params.ct} --id ${params.samplename} --outdir ${params.outdir}
    kallisto index -i ${params.transindex} -k 31 --make-unique ${transcriptome_file}
    """
}

/*
 * Step 2. Do pseudoalignment with kallisto for plate-based data
 */
process pseudoalPlate {
    storeDir params.outdir

    input:
    file index from transcriptome_index
    file batchkal from batch_kal.collect()

    output:
    file params.samplename into kallisto_pseudo

    when: params.protocol=='plate'

    script:
    """
    kallisto pseudo -t ${params.cores} --quant \\
        -i $index \\
        -o ${params.samplename} \\
        -b $batchkal
    """
}

/*
 * Step 2. Do pseudoalignment with kallisto
 */
process pseudoal {

    //publishDir params.outdir, mode: 'copy'
    storeDir params.outdir

    input:
    file index from transcriptome_index
    file reads from read_files_kallisto.collect()

    output:
    file params.samplename into kallisto_bus_to_sort
    file "${params.samplename}/output.bus"
    file "${params.samplename}/matrix.ec"
    file "${params.samplename}/transcripts.txt"

    when: params.protocol!='plate'

    script:
    if(params.protocol=='visiumv1')
        """
        kallisto bus \\
            -i $index \\
            -o ${params.samplename} \\
            -x 0,0,16:0,16,28:1,0,0 \\
            -t ${params.cores} \\
            $reads
        """
    else if(params.protocol=='sc5pe')
        """
        kallisto bus \\
            -i $index \\
            -o ${params.samplename} \\
            -x 0,0,16:0,16,26:0,26,0,1,0,0 \\
            -t ${params.cores} \\
            $reads
        """
    else if(params.protocol=='sc5pe')
        """
        kallisto bus \\
            -i $index \\
            -o ${params.samplename} \\
            -x 0,0,16:0,16,26:0,26,0,1,0,0 \\
            -t 16 \\
            $reads
        """
    else
        """
        kallisto bus \\
            -i $index \\
            -o ${params.samplename} \\
            -x ${params.protocol} \\
            -t ${params.cores} \\
            $reads
        """
}
