/*
 * Example commands:
 *
 * nextflow run kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/10xv3_whitelist.txt --samplename "test_10xv3" --outdir ./ --protocol 10xv3 --reads "test_datasets/10xv3/Brain_Tumor_3p_fastqs/*.fastq.gz" --cores 4
 *
 * nextflow run kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/barcode_whitelists/visium-v1_whitelist_kallisto.txt --samplename "test_visium_FFPE" --outdir ./ --protocol visiumv1 --reads "test_datasets/visium_FFPE/Visium_FFPE_Human_Ovarian_Cancer_fastqs/*.fastq.gz" --images "V10L13-020" --imagear "D1" --imagef test_datasets/visium_FFPE/Visium_FFPE_Human_Ovarian_Cancer_image.jpeg --cores 2
 */

// sc5pe definition taken from https://github.com/pachterlab/kallisto/issues/287

/*
 * Defines necessary parameters
 * reads is used for path for fastq files, but also for the batch file for plate-based and bulk data
 * reads must be in the format *R1.* and *R2.*
 */
// output parameters
params.samplename = false
params.outdir = "./"
// transcriptome index parameters
params.transcriptome = false
params.transindex = false
// processing parameters
params.protocol = false 
params.reads = false // all reads or batch file
params.white = false 
params.t2g = false
params.cores = 1
// Visium image parameters
params.imageal = false
params.imagef = false
params.imagear = false
params.images = false
// RNA velocity parameters
params.velomode = false
params.genome = false
params.gtf = false
params.overhang = 0 // intron overhangs should be the length of the biological read minus 1


if(false){
/*
 * Input parameters validation and preparation
 */
//----- TRANSCRIPTOME -----//
//if(!params.transindex) exit 1, "Transcriptome index file path is mandatory (will be created if it does not exist)."
transindex_file = file(params.transindex)
if(!transindex_file.exists() and params.transcriptome!="") transcriptome_file = file(params.transcriptome)
else exit 1, "Missing transcriptome file or transcriptome index file."

//----- READS -----//
R1files = Channel
    .from(params.reads.tokenize())
    .flatMap{ files(it) }
    .filter{it =~ /.*R1.*/}
    .toSortedList()
    .flatten()
    //.view()
R2files = Channel
    .from(params.reads.tokenize())
    .flatMap{ files(it) }
    .filter{it =~ /.*R2.*/}
    .toSortedList()
    .flatten()
    //.view()
R1files
    .merge(R2files){ a, b -> tuple(a,b) }
    .flatten()
    //.view()
    .set{read_files_kallisto}
if(params.protocol=='batch'){
    Channel.fromPath(params.reads)
        .set{batch_kal}
} else{ batch_kal = "" }

//----- OTHER FILES -----//
if(!params.samplename) exit 1, "Please provide a name for this sample"

if(!params.protocol) exit 1, "Please provide an adequate protocol"

if(!params.white && (params.protocol!='batch' && params.protocol!='bulk_quant')){
    exit 1, "Barcode whitelist is mandatory for Chromium, Visium, and ParseBio runs."
}

if(!params.t2g){
    exit 1, "Transcriptome-to-gene reference is required for quantification."
} else{
    Channel.fromPath(params.t2g)
        .set{t2g_kal}
}
}
//----- RNA VELOCITY -----//
gtffname = file(params.gtf)["baseName"].replaceFirst(/.gtf/, "")



/*
 * For some reason pseudoalignment doesn't work if it is the first step
 */
process starting{

    script:
    """
    echo Pseudoalignment doesn't work as a first step, so we'll have this here.
    """
}

/*
 * Builds the transcriptome index, if it doesn't exist
 */
process index {

    input:
    path transcriptome_file

    output:
    path("${params.transindex}")
    
    storeDir transcriptome_file.getParent()

    script:
    """
    echo Transcriptome index will be stored in the parent directory of the reference.
    kallisto index -i ${params.transindex} -k 31 -t ${params.cores} --make-unique ${transcriptome_file}
    """
}

/*
 * Do pseudoalignment with kallisto for batch data
 */
process pseudoalBatch {

    input:
    path index
    path batchkal

    output:
    path "${params.samplename}/output.bus"
    path "${params.samplename}/matrix.ec"
    path "${params.samplename}/transcripts.txt"
    path "${params.samplename}/run_info.json"
    
    storeDir params.outdir

    when: params.protocol=='batch'

    script:
    """
    kallisto bus -t ${params.cores} -i $index -o ${params.samplename} \\
        --batch-barcodes -B $batchkal
    """
}

/*
 * Do pseudoalignment with kallisto for bulk data, one dataset
 */
process bulk_quant {

    input:
    path index
    path reads

    output:
    path "${params.samplename}"
    
    storeDir params.outdir

    when: params.protocol=='bulk_quant'

    script:
    """
    kallisto quant -t ${params.cores} -i $index -o ${params.samplename} ${reads}
    """
}

/*
 * Do pseudoalignment with kallisto
 */
process pseudoal {
  
    input:
    path index
    path reads

    output:
    path "${params.samplename}/output.bus"
    path "${params.samplename}/matrix.ec"
    path "${params.samplename}/transcripts.txt"
    path "${params.samplename}/run_info.json"
    
    storeDir params.outdir
    
    when: params.protocol!='batch' && params.protocol!='bulk_quant'

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

/*
 * Correct barcodes and sort
 * A future version of kallisto may include options to correct the UMI
 */
process corrsort {

    input:
    path outbus
    path white

    output:
    path("${params.samplename}/output.cor.sort.bus")
    
    storeDir params.outdir

    when: params.protocol!='bulk_quant'

    script:
    if(params.protocol=='batch')
    """
    mkdir -p ${params.samplename}
    bustools sort -o ${params.samplename}/output.cor.sort.bus -t ${params.cores} ${outbus}
    """
    else
    """
    mkdir -p ${params.samplename}
    bustools correct -w ${white} -o ${params.samplename}/output.cor.bus ${outbus}
    bustools sort -o ${params.samplename}/output.cor.sort.bus \\
        -t ${params.cores} ${params.samplename}/output.cor.bus
    """
}

/*
 * Obtain the UMI read counts
 */
process umicounts {

    input:
    path outbus

    output:
    path "${params.samplename}/umicount.txt.gz"
    
    storeDir params.outdir

    when: params.protocol!='batch' && params.protocol!='bulk_quant'

    script:
    """
    mkdir -p ${params.samplename}
    bustools text -o ${params.samplename}/umicount.txt ${outbus}
    gzip ${params.samplename}/umicount.txt
    """
}

/*
 * Obtain the counts
 */
process countbus {
  
    input:
    path outbus
    path outmat
    path outtrans
    path t2g

    output:
    path "${params.samplename}/genecounts.mtx"
    path "${params.samplename}/genecounts.barcodes.txt"
    path "${params.samplename}/genecounts.genes.txt"

    storeDir params.outdir
    
    when: params.protocol!='bulk_quant'

    script:
    """
    mkdir -p ${params.samplename}
    bustools count --em -t ${outtrans} -e ${outmat} -g $t2g --genecounts \\
    -o ${params.samplename}/genecounts ${outbus}
    """
}

/*
 * Make Seurat object for plate
 */
process makeSeuratPlate {

    input:
    path outmtx
    path outbc
    path outg

    output:
    path "${params.samplename}_srat.RDS"
    
    storeDir "${params.outdir}/${params.samplename}"

    when: params.protocol=='batch' && params.protocol!='bulk_quant'

    script:
    """
    #!Rscript --vanilla

    library(Seurat)

    # read in data
    topdir = "${params.samplename}" # source dir
    exp = Matrix::readMM("${outmtx}") #read matrix
    bc = read.csv("${outbc}", header = F, stringsAsFactors = F)
    g = read.csv("${outg}", header = F, stringsAsFactors = F)
    dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
    count.data = Matrix::t(exp)
    rm(exp) # help with memory management
    gc()

    # create Seurat object
    srat = CreateSeuratObject(counts = count.data)

    # get MT% (genes curated from NCBI chrMT genes)
    mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
                "ATP8", "MT-CO1", "COI", "LOC9829747")
    mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
    mtall = paste0("-",paste(mtgenes, collapse="\$|-"))
    mtgenes = grep(mtall, rownames(srat), value = T, ignore.case = T)
    #mtgenes = mtgenes[mtgenes %in% g[,1]]
    srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "RNA",
                                features = mtgenes)

    saveRDS(srat, file = "${params.samplename}_srat.RDS")
    """

}

/*
 * Make Seurat object for 10x
 */
process makeSeurat10x {

    input:
    path outmtx
    path outbc
    path outg
    path umic

    output:
    path "${params.samplename}_UMIrank.pdf"
    path "${params.samplename}_UMIduplication.pdf"
    path "${params.samplename}_srat.RDS"
    
    storeDir "${params.outdir}/${params.samplename}"

    when: params.protocol=='10xv3' || params.protocol=='10xv2' || params.protocol=='sc5pe'

    script:
    """
    #!Rscript --vanilla

    library(Seurat)
    library(DropletUtils)
    library(data.table)

    # read in data
    topdir = "${params.samplename}" # source dir
    exp = Matrix::readMM("${outmtx}") #read matrix
    bc = read.csv("${outbc}", header = F, stringsAsFactors = F)
    g = read.csv("${outg}", header = F, stringsAsFactors = F)
    dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
    count.data = Matrix::t(exp)
    rm(exp) # help with memory management
    gc()

    # get emptyDrops and default cutoff cell estimates
    iscell_dd = defaultDrops(count.data, expected = 5000)
    eout = emptyDrops(count.data, lower = 200)
    eout\$FDR[is.na(eout\$FDR)] = 1
    iscell_ed = eout\$FDR<=0.01
    meta = data.frame(row.names = paste0(bc\$V1,"-1"),
                      iscell_dd = iscell_dd, iscell_ed = iscell_ed)

    # plot rankings for number of UMI
    br.out <- barcodeRanks(count.data)
    pdf("${params.samplename}_UMIrank.pdf", height = 5, width = 5, useDingbats = F)
    plot(br.out\$rank, br.out\$total, log="xy", xlab="Rank", ylab="Total")
    o <- order(br.out\$rank)
    lines(br.out\$rank[o], br.out\$fitted[o], col="red")
    abline(h=metadata(br.out)\$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)\$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
        legend=c("knee", "inflection"))
    dev.off()

    # UMI duplication
    umi = fread("${umic}", sep = "\t", header = F, stringsAsFactors = F)
    sumUMI = c()
    sumi = sum(umi\$V4)
    for(i in 0:250){ sumUMI = c(sumUMI, sum(umi\$V4[umi\$V4>i])/sumi) }
    pdf("${params.samplename}_UMIduplication.pdf", height = 3.5, width = 7, useDingbats = F)
    par(mfrow = c(1,2))
    plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
         xlab = "More than xx UMI", main = "${params.samplename}")
    diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
    plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
         xlab = "More than xx UMI", main = "${params.samplename}")
    dev.off()
    rm(umi) # help with memory management
    gc()

    # create Seurat object
    ## we're only keeping what might potentially be a cell (by DD or ED)
    srat = CreateSeuratObject(counts = count.data[,iscell_dd | iscell_ed],
                              meta.data = meta[iscell_dd | iscell_ed,])
    amb_prop = estimateAmbience(count.data)[rownames(srat@assays\$RNA@meta.features)]
    srat@assays\$RNA@meta.features = data.frame(row.names = rownames(srat@assays\$RNA@meta.features),
                                                "ambient_prop" = amb_prop)

    # get MT% (genes curated from NCBI chrMT genes)
    mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
                "ATP8", "MT-CO1", "COI", "LOC9829747")
    mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
    mtall = paste0("-",paste(mtgenes, collapse="\$|-"))
    mtgenes = grep(mtall, rownames(srat), value = T, ignore.case = T)
    #mtgenes = mtgenes[mtgenes %in% g[,1]]
    srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "RNA",
                                features = mtgenes)

    saveRDS(srat, file = "${params.samplename}_srat.RDS")
    """

}

/*
 * Make Seurat object for ParseBio
 */
process makeSeuratParse {

    input:
    path outmtx
    path outbc
    path outg
    path umic

    output:
    path "${params.samplename}_UMIrank.pdf"
    path "${params.samplename}_srat.RDS"
    
    storeDir "${params.outdir}/${params.samplename}"

    when: params.protocol=='SPLiT-seq'

    script:
    """
    #!Rscript --vanilla

    library(Seurat)
    library(DropletUtils)
    library(DelayedArray)
    library(data.table)

    # read in data
    topdir = "${params.samplename}" # source dir
    exp = Matrix::readMM("${outmtx}") #read matrix
    bc = read.csv(paste0(topdir, "/genecounts.barcodes.txt"), header = F, stringsAsFactors = F)
    g = read.csv(paste0(topdir, "/genecounts.genes.txt"), header = F, stringsAsFactors = F)
    dimnames(exp) = list(paste0(bc[,1],"-1"), g[,1]) # number added because of seurat format for barcodes
    exp = exp[Matrix::rowSums(exp)>0,]

    # summarise cells based on well
    well_file = "${params.white}"
    well_file = gsub("_whitelist_", "_wells_", well_file)
    well_file = read.csv(well_file, header = T)
    well_file\$barcode = paste0(well_file\$barcode, "-1")

    maxbc1 = max(well_file\$bci1)/2
    l_sum = list()
    for(i in 1:maxbc1){
      sub_well_file = well_file[well_file\$bci1 %in% c(i, i+maxbc1),]
      sub_well_file = sub_well_file[sub_well_file\$barcode %in% rownames(exp),]

      sub_exp = as(exp[sub_well_file\$barcode,], "CsparseMatrix")
      sub_exp = DelayedArray::rowsum(sub_exp, group = sub_well_file\$well_combo[match(rownames(sub_exp), sub_well_file\$barcode)])
      rownames(sub_exp) = sub_well_file\$barcode[!duplicated(sub_well_file[,2])]

      l_sum[[i]] = as(sub_exp, "CsparseMatrix")
    }
    exp = do.call(rbind, l_sum)
    sub_well_file = well_file[well_file\$barcode %in% rownames(exp),]

    # transpose
    count.data = Matrix::t(exp)

    # plot rankings for number of UMI
    br.out <- barcodeRanks(count.data)
    pdf("${params.samplename}_UMIrank.pdf", height = 5, width = 5, useDingbats = F)
    plot(br.out\$rank, br.out\$total, log="xy", xlab="Rank", ylab="Total")
    o <- order(br.out\$rank)
    lines(br.out\$rank[o], br.out\$fitted[o], col="red")
    abline(h=metadata(br.out)\$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)\$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
        legend=c("knee", "inflection"))
    dev.off()

    # get emptyDrops and default cutoff cell estimates
    iscell_dd = defaultDrops(count.data, expected = 10000)
    eout = emptyDrops(count.data, lower = 200)
    eout\$FDR[is.na(eout\$FDR)] = 1
    iscell_ed = eout\$FDR<=0.01
    meta = data.frame(row.names = colnames(count.data),
                      iscell_dd = iscell_dd,
                      iscell_ed = iscell_ed,
                      iscell_inflexion = br.out\$total>=metadata(br.out)\$inflection,
                      iscell_knee = br.out\$total>=metadata(br.out)\$knee)
    meta = merge(meta, sub_well_file, by.x = 0, by.y = 1)
    rownames(meta) = meta[,1]
    meta = meta[,-1]

    # UMI duplication
    #umi = fread("${umic}", sep = "\t", header = F, stringsAsFactors = F)
    #sumUMI = c()
    #sumi = sum(umi\$V4)
    #for(i in 0:250){ sumUMI = c(sumUMI, sum(umi\$V4[umi\$V4>i])/sumi) }
    #pdf("${params.samplename}_UMIduplication.pdf", height = 3.5, width = 7, useDingbats = F)
    #par(mfrow = c(1,2))
    #plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
    #     xlab = "More than xx UMI", main = "${params.samplename}")
    #diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
    #plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
    #     xlab = "More than xx UMI", main = "${params.samplename}")
    #dev.off()
    #rm(umi) # help with memory management
    #gc()

    # create Seurat object
    ## we're only keeping what might potentially be a cell (by DD or ED)
    cells_keep = meta\$iscell_dd | meta\$iscell_ed | meta\$iscell_inflexion | meta\$iscell_knee
    sub_meta = meta[cells_keep,]
    srat = CreateSeuratObject(counts = count.data[,rownames(sub_meta)],
                              meta.data = sub_meta)
    amb_prop = estimateAmbience(count.data)[rownames(srat@assays\$RNA@meta.features)]
    srat@assays\$RNA@meta.features = data.frame(row.names = rownames(srat@assays\$RNA@meta.features),
                                                "ambient_prop" = amb_prop)

    # get MT% (genes curated from NCBI chrMT genes)
    mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
                "ATP8", "MT-CO1", "COI", "LOC9829747")
    mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
    mtall = paste0("-",paste(mtgenes, collapse="\$|-"))
    mtgenes = grep(mtall, rownames(srat), value = T, ignore.case = T)
    #mtgenes = mtgenes[mtgenes %in% g[,1]]
    srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "RNA",
                                features = mtgenes)

    saveRDS(srat, file = "${params.samplename}_srat.RDS")
    """

}

/*
 * Get tissue alignment and fiducials
 */
process getTissue {
    
    input:
    val imageal
    val imagef
    val imagear
    val images

    output:
    path "${params.samplename}/outs/spatial"
    //path "${params.samplename}/outs/spatial/aligned_fiducials.jpg"
    //path "${params.samplename}/outs/spatial/scalefactors_json.json"
    //path "${params.samplename}/outs/spatial/tissue_lowres_image.png"
    //path "${params.samplename}/outs/spatial/detected_tissue_image.jpg"
    //path "${params.samplename}/outs/spatial/tissue_hires_image.png"
    //path "${params.samplename}/outs/spatial/tissue_positions.csv"
    
    storeDir "${params.outdir}/${params.samplename}"

    when: params.protocol=='visiumv1'

    script:
    if(!params.imageal)
      """
      spaceranger count --id=${params.samplename} --fastqs=/opt/mock_fastq \\
      --transcriptome=/opt/eGFP --image=${imagef} --slide=${images} \\
      --area=${imagear} --localcores=${params.cores}
      """
    else
      """
      spaceranger count --id=${params.samplename} --fastqs=/opt/mock_fastq \\
      --transcriptome=/opt/eGFP --image=${imagef} --slide=${images} \\
      --area=${imagear} --loupe-alignment=${imageal} --localcores=${params.cores}
      """

}

/*
 * Make Seurat object for Visium v1
 */
process makeSeuratVisium {

    input:
    path outmtx
    path outbc
    path outg
    path umic
    // from getTissue
    path spatial
    //path fiducials
    //path scalefactors
    //path tissue_lowres
    //path detected_tissue
    //path tissue_hires
    //path tissue_positions

    output:
    path "${params.samplename}_UMIrank.pdf"
    path "${params.samplename}_UMIduplication.pdf"
    path "${params.samplename}_srat.RDS"
    
    storeDir "${params.outdir}/${params.samplename}"

    when: params.protocol=='visiumv1'

    script:
    """
    #!Rscript --vanilla

    library(Seurat)
    library(DropletUtils)
    library(data.table)

    # read in data
    topdir = "${params.samplename}" # source dir
    exp = Matrix::readMM("${outmtx}") #read matrix
    bc = read.csv("${outbc}", header = F, stringsAsFactors = F)
    g = read.csv("${outg}", header = F, stringsAsFactors = F)
    dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
    count.data = Matrix::t(exp)
    rm(exp) # help with memory management
    gc()

    # plot rankings for number of UMI
    br.out <- barcodeRanks(count.data)
    pdf("${params.samplename}_UMIrank.pdf", height = 5, width = 5, useDingbats = F)
    plot(br.out\$rank, br.out\$total, log="xy", xlab="Rank", ylab="Total")
    o <- order(br.out\$rank)
    lines(br.out\$rank[o], br.out\$fitted[o], col="red")
    abline(h=metadata(br.out)\$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)\$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
        legend=c("knee", "inflection"))
    dev.off()

    # UMI duplication
    umi = fread("${umic}", sep = "\t", header = F, stringsAsFactors = F)
    sumUMI = c()
    sumi = sum(umi\$V4)
    for(i in 0:250){ sumUMI = c(sumUMI, sum(umi\$V4[umi\$V4>i])/sumi) }
    pdf("${params.samplename}_UMIduplication.pdf", height = 3.5, width = 7, useDingbats = F)
    par(mfrow = c(1,2))
    plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
         xlab = "More than xx UMI", main = "${params.samplename}")
    diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
    plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
         xlab = "More than xx UMI", main = "${params.samplename}")
    dev.off()
    rm(umi) # help with memory management
    gc()

    # create Seurat object
    srat <- CreateSeuratObject(counts = count.data, assay = "Spatial",
                               meta.data = read.csv("${spatial}/tissue_positions.csv", 
                                                    header = T, row.names = 1)) # create object
    image <- Read10X_Image(image.dir = "${spatial}",
                           filter.matrix = FALSE) # read in the images
    image <- image[Cells(x = srat)] # filter image by the spots
    DefaultAssay(object = image) <- "Spatial" # set default assay
    srat[["slice1"]] <- image # slice name might be changed

    # estimate ambient proportion
    amb_prop = estimateAmbience(count.data)[rownames(srat@assays\$Spatial@meta.features)]
    srat@assays\$Spatial@meta.features = data.frame(row.names = rownames(srat@assays\$Spatial@meta.features),
                                                    "ambient_prop" = amb_prop)

    # get MT% (genes curated from NCBI chrMT genes)
    mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
                "ATP8", "MT-CO1", "COI", "LOC9829747")
    mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
    mtall = paste0("-",paste(mtgenes, collapse="\$|-"))
    mtgenes = grep(mtall, rownames(srat), value = T, ignore.case = T)
    #mtgenes = mtgenes[mtgenes %in% g[,1]]
    srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "Spatial",
                                features = mtgenes)

    saveRDS(srat, file = "${params.samplename}_srat.RDS")
    """
}


/*
 * RNA VELOCITY PROCESSES
 */
// parameters (defined above)
// params.velomode = false
// params.genome = false
// params.gtf = false
// params.overhang = 0
/*
 * Extract Introns and Exons from genome and GTF
 */
process getIntronCoord {
  
    input:
    path gtf
    val overhang
    
    output:
    path "${gtffname}_introns.gtf"
    path "${gtffname}_intronsPlus${overhang}.gtf"
    
    storeDir file(params.gtf).getParent()
    
    when: params.velomode
    
    script:
    """
    #!Rscript --vanilla
    
    library(gread)
    library(rtracklayer)
    
    gtf = read_format("${gtf}")
    ans = construct_introns(gtf, update=F)[]
    ans = sort(ans)
    ans\$transcript_id = paste0(ans\$transcript_id, ".I", 1:length(ans))
    
    ans_plus = ans+${overhang}
    start(ans_plus)[start(ans_plus)<1] = 1
    
    export(ans, paste0("${gtffname}", "_introns.gtf"))
    export(ans_plus, paste0("${gtffname}", "_intronsPlus", ${overhang}, ".gtf"))
    """
}

/*
 * Make a reference consisting of cDNA transcripts and introns
 */
process intronTranscriptRef {
    input:
    path gtf
    path genome
    path transcriptome
    
    output:
    path "${gtffname}_intronsPlus${params.overhang}_corr.gtf"
    path "${genome}.fai"
    path "${gtffname}_intronsPlus${params.overhang}_corr.fa"
    path "transcriptome_and_intronsPlus${params.overhang}_corr.fa"
    
    storeDir file(params.genome).getParent()
    
    when: params.velomode
    
    script:
    """
    sed 's/sequence_feature/exon/g' ${gtf} | sed 's/intron/exon/g' - > ${gtffname}_intronsPlus${params.overhang}_corr.gtf
    samtools faidx ${genome}
    gffread -g ${genome} -w ${gtffname}_intronsPlus${params.overhang}_corr.fa ${gtffname}_intronsPlus${params.overhang}_corr.gtf
    gzip -cd ${transcriptome} | cat - ${gtffname}_intronsPlus${params.overhang}_corr.fa > transcriptome_and_intronsPlus${params.overhang}_corr.fa
    """
}

/*
 * Builds the transcriptome index for introns and exons
 */
process indexInEx {

    input:
    path transcriptome_introns

    output:
    path "transcriptome_and_intronsPlus${params.overhang}_corr.kalid"
    
    storeDir file(params.genome).getParent()

    script:
    """
    kallisto index -i transcriptome_and_intronsPlus${params.overhang}_corr.kalid \\
        -k 31 -t ${params.cores} --make-unique ${transcriptome_introns}
    """
}

/*
 * Subset introns and exons
 */
process captureInEx {

    input:
    path outbus // corrsort.out
    path matrix // pseudoal.out[1]
    path transcripts // pseudoal.out[2]
    path intronsfile
    path exonsfile
    // don't forget this works with complement sets!
    // https://bustools.github.io/BUS_notebooks_R/velocity.html

    output:
    path "${params.samplename}/spliced.bus"
    path "${params.samplename}/unspliced.bus"
    
    storeDir params.outdir
    
    when: params.velomode

    script:
    """
    bustools capture -s -x -o ${params.samplename}/spliced.bus -c $intronsfile \\
        -e ${matrix} -t ${transcripts} ${outbus}
    bustools capture -s -x -o ${params.samplename}/unspliced.bus -c $exonsfile \\
        -e ${matrix} -t ${transcripts} ${outbus}
    """
}

/*
 * Count introns and exons
 */
process countInEx {
    input:
    
    
    output:
        
    storeDir params.outdir
    
    when: params.velomode
    
    script:
    """
    
    """
}

/*
 * Process introns and exons
 */
process processInEx {
    input:
    
    
    output:
        
    storeDir params.outdir
    
    when: params.velomode

    script:
    """
    
    """
}



/*
 * Execute workflow
 */
workflow {
  starting()
  if(!params.velomode){
    // indexing
    index(file(params.transcriptome))
    if(params.protocol=='bulk_quant'){
      bulk_quant(index.out, read_files_kallisto.collect())
    } else{
      // pseudoalignment
      pseudoalBatch(index.out, batch_kal.collect())
      pseudoal(index.out, read_files_kallisto.collect())
      // sorting and barcode corrections
      if(params.protocol=='batch'){
        corrsort(pseudoalBatch.out[0], pseudoalBatch.out[1])
        // pseudoalBatch.out[1] only serves as placeholder
      } else{
        corrsort(pseudoal.out[0], file(params.white))
      }
      // umi counts
      umicounts(corrsort.out)
      // counts per gene
      if(params.protocol=='batch'){
        countbus(corrsort.out, pseudoalBatch.out[1], pseudoalBatch.out[2], t2g_kal.collect())
      } else{
        countbus(corrsort.out, pseudoal.out[1], pseudoal.out[2], t2g_kal.collect())
      }
      // tissue image
      getTissue(params.imageal, params.imagef, params.imagear, params.images)
      // make Seurat
      makeSeuratPlate(countbus.out[0], countbus.out[1], countbus.out[2])
      makeSeurat10x(countbus.out[0], countbus.out[1], countbus.out[2], umicounts.out)
      makeSeuratParse(countbus.out[0], countbus.out[1], countbus.out[2], umicounts.out)
      makeSeuratVisium(countbus.out[0], countbus.out[1], countbus.out[2], umicounts.out,
      getTissue.out[0])
    }
  } else{
    getIntronCoord(file(params.gtf), params.overhang)
    intronTranscriptRef(getIntronCoord.out[1], file(params.genome), file(params.transcriptome))
    indexInEx(intronTranscriptRef.out[3])
  }
}


