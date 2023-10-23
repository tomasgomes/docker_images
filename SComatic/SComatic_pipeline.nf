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
params.maxcov = 150
params.cores = 12
params.cellsites = true
params.cellgenotypes = true



/*
 * Input parameters validation and preparation
 */
if(!params.bam) exit 1, "Please provide a sorted BAM file."
if(!params.ct) exit 1, "Please provide a file with cell type annotation. Should be a tab-separated table with Index and Cell_type columns."
if(!params.samplename) exit 1, "Please provide a name for this sample."
if(!params.genome) exit 1, "Please provide a genome reference in FASTA format."



/*
 * Step 0. Get things going
 */
process starting {

    script:
    """
    ls -lh
    """
}

/*
 * Step 1. Split the BAM file by cell type
 */
process splitBAM {

  input:
  val inbam // has to be val because of the index file
  path metact

  output:
  path("${params.samplename}.report.txt")
  path("${params.samplename}.*.bam")
  path("${params.samplename}.*.bam.bai")
  
  storeDir "${params.outdir}"

  script:
  """
  mkdir -p /rootvol/${params.outdir}
  python3 /rootvol/SComatic/scripts/SplitBam/SplitBamCellTypes.py --bam ${inbam} \\
      --meta ${metact} --id ${params.samplename} --outdir ./
  """
}

/*
 * Step 2. Detect base variations per cell
 */
process baseCellCounts {
  
  input:
  path subbam // has to be val because of the index file
  path subbai

  output:
  path outtsv

  storeDir "${params.outdir}/tables"
  
  script:
  outtsv = subbam.name.replaceAll("bam","tsv")
  """
  mkdir -p ${params.outdir}/tables
  python3 /rootvol/SComatic/scripts/BaseCellCounter/BaseCellCounter.py \\
      --bam ${subbam} --ref ${params.genome} --chrom all --out_folder ./ \\
      --nprocs ${params.cores} --tmp_dir ./${params.samplename}_tmp
  """
}

/*
 * Step 3. Merge base counts
 */
process mergeCounts {
  
  input:
  path tables_calls
  
  output:
  path("merged_tables.tsv")
  
  storeDir "${params.outdir}/tables"
  
  script:
  """
  python3 /rootvol/SComatic/scripts/MergeCounts/MergeBaseCellCounts.py \\
      --tsv_folder ./ --outfile ./merged_tables.tsv
  """
}

/*
 * Step 4. Call variation per cell type
 */
process baseCallStep1 {
  
  input:
  path merged_tables
  
  output:
  path "basecalling1.calling.step1.tsv"
  
  storeDir "${params.outdir}"
  
  script:
  """
  python3 /rootvol/SComatic/scripts/BaseCellCalling/BaseCellCalling.step1.py \\
      --infile ${merged_tables} --outfile ./basecalling1 \\
      --ref ${params.genome} --min_cells ${params.mincells}
  """
}

/*
 * Step 5. Compare calls to references
 */
process baseCallStep2 {
  
  input:
  path basecalling1
  
  output:
  path "basecalling2.calling.step2.tsv"
  
  storeDir "${params.outdir}"
  
  script:
  """
  python3 /rootvol/SComatic/scripts/BaseCellCalling/BaseCellCalling.step2.py \\
      --infile ${basecalling1} --outfile basecalling2 \\
      --editing /rootvol/SComatic/RNAediting/AllEditingSites.hg38.txt \\
      --pon /rootvol/SComatic/PoNs/PoN.scRNAseq.hg38.tsv
  """
}

/*
 * Step 6. Collect all variation sites
 */
process callableSites {
  
  input:
  file basecalling2
  
  output:
  file "callableSites.coverage_cell_count.report.tsv"
  
  storeDir "${params.outdir}"
  
  script:
  """
  python3 /rootvol/SComatic/scripts/GetCallableSites/GetAllCallableSites.py \\
      --infile ${basecalling2} --outfile callableSites --max_cov ${params.maxcov} \\
      --min_cell_types ${params.mincelltypes}
  """
}

/*
 * Step 7. Determine the callable sites per cell
 */
process uniqueCallableSites {
  
  input:
  path basecalling1
  path subbam
  path subbai
  
  output:
  path outtsv2

  storeDir "${params.outdir}/UniqueCellCallableSites"
  
  when: params.cellsites==true
  
  script:
  outtsv2 = subbam.name.replaceAll("bam","SitesPerCell.tsv")
  ctn = subbam.name.tokenize(".")[1]
  """
  mkdir -p ${params.outdir}/UniqueCellCallableSites
  mkdir -p ./${params.samplename}_${ctn}_tmp
  python3 /rootvol/SComatic/scripts/SitesPerCell/SitesPerCell.py \\
      --bam ${subbam} --infile ${basecalling1} --ref ${params.genome} \\
      --out_folder ./ --tmp_dir ./${params.samplename}_${ctn}_tmp --nprocs ${params.cores}
  """
}

/*
 * Step 8. Get genotypes for each cell
 */
process cellGenotypes {
  
  input:
  path basecalling2
  path subbam
  path subbai
  
  output:
  path("${params.samplename}.${ctn}.single_cell_genotype.tsv")

  storeDir "${params.outdir}/SingleCellAlleles"
  
  when: params.cellgenotypes==true
  
  script:
  ctn = subbam.name.tokenize(".")[1]
  """
  mkdir -p ${params.outdir}/SingleCellAlleles
  mkdir -p ./${params.samplename}_${ctn}_tmp
  python3 /rootvol/SComatic/scripts/SingleCellGenotype/SingleCellGenotype.py \\
      --bam ${subbam} --infile ${basecalling2} --nprocs ${params.cores} \\
      --meta ${params.ct} --ref ${params.genome} --tmp_dir ./${params.samplename}_${ctn}_tmp \\
      --outfile ${params.samplename}.${ctn}.single_cell_genotype.tsv \\
  """
}



/*
 * Execute workflow
 */
workflow {
  starting()
  splitBAM(params.bam, params.ct)
  splitBAM.out[1].view()
  baseCellCounts(splitBAM.out[1] | flatten, splitBAM.out[2] | flatten)
  baseCellCounts.out.collect().view()
  mergeCounts(baseCellCounts.out.collect())
  mergeCounts.out.view()
  baseCallStep1(mergeCounts.out)
  baseCallStep2(baseCallStep1.out)
  callableSites(baseCallStep2.out)
  uniqueCallableSites(baseCallStep1.out, splitBAM.out[1] | flatten, splitBAM.out[2] | flatten)
  cellGenotypes(baseCallStep2.out, splitBAM.out[1] | flatten, splitBAM.out[2] | flatten)
}





