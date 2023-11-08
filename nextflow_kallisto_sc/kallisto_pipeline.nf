/*
 * Example command:
 *
 * nextflow run kallisto_pipeline.nf --transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial.fa --transindex AmexT_v47_artificial.kalid --t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial_genenames_t2g.txt --white /links/groups/treutlein/USERS/tomasgomes/gene_refs/other/10xv2_whitelist.txt --samplename "GER006_10x" --outdir /links/groups/treutlein/USERS/tomasgomes/data/axolotl/ --protocol 10xv2 --reads "/links/groups/treutlein/USERS/tomasgomes/projects/axolotl/data/raw/Gerber_all10x/GER006_10x/*.fastq.gz"
 *
 * nextflow run kallisto_pipeline.nf --transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial.fa --transindex AmexT_v47_artificial.kalid --t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial_genenames_t2g.txt --samplename "Gerber_plate" --outdir /links/groups/treutlein/USERS/tomasgomes/data/axolotl/ --protocol plate --reads /links/groups/treutlein/USERS/tomasgomes/projects/axolotl/data/raw/Gerber_allcells/kallisto_batch.txt
 */

// sc5pe taken from https://github.com/pachterlab/kallisto/issues/287

/*
 * Defines necessary parameters
 * reads is used for path for fastq files, but also for the batch file for plate-based data
 * reads must be in the format *R1* and *R2*
 */
params.samplename = false //params.samplename = "GER006_10x"
params.outdir = "./"
params.protocol = false //params.protocol = "10xv2"
params.transcriptome = false //params.transcriptome = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial.fa"
params.transindex = false //params.transindex = "AmexT_v47_artificial.kalid"
params.reads = false //params.reads = "/links/groups/treutlein/USERS/tomasgomes/projects/axolotl/data/raw/Gerber_all10x/GER006_10x/*.fastq.gz"
params.white = false //params.white = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/other/10xv2_whitelist.txt"
params.t2g = false //params.t2g = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial_genenames_t2g.txt"
params.imageal = false
params.imagef = false
params.cores = 12



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
if(params.protocol=='plate'){
    Channel.fromPath(params.reads)
        .set{batch_kal}
} else{ batch_kal = "" }

//----- OTHER FILES -----//
if(!params.samplename) exit 1, "Please provide a name for this sample"

if(!params.protocol) exit 1, "Please provide an adequate protocol"

if(!params.white && params.protocol!='plate'){
    exit 1, "Barcode whitelist is mandatory for Chromium, Visium, and ParseBio runs."
} else if (params.protocol!='plate'){
    Channel.fromPath(params.white)
        .set{bc_wl_kal}
} else{ bc_wl_kal = "" }

if(!params.t2g && params.protocol!='plate'){
    exit 1, "Transcriptome-to-gene reference is required for quantification."
} else if(params.protocol!='plate'){
    Channel.fromPath(params.t2g)
        .set{t2g_kal}
    t2g_plate = ""
} else{
    Channel.fromPath(params.t2g)
        .set{t2g_plate}
    t2g_kal = ""
}

if(file("${params.outdir}${params.samplename}").exists()) println "Output folder already exists. No files will be overwritten, but execution may fail."


/*
 * Step 0. For some reason pseudoalignment doesn't work if it is the first step
 */
process starting{

    script:
    """
    echo Pseudoalignment doesn't work as a first step, so we'll have this here.
    """
}

/*
 * Step 1. Builds the transcriptome index, if it doesn't exist
 */
process index {

    input:
    val transcriptome_file

    output:
    path("${params.transindex}")
    
    storeDir transcriptome_file.getParent()

    script:
    """
    echo Transcriptome index will be stored in the parent directory of the reference.
    kallisto index -i ${params.transindex} -k 31 --make-unique ${transcriptome_file}
    """
}

/*
 * Step 2. Do pseudoalignment with kallisto for plate-based data
 */
process pseudoalPlate {

    input:
    path index
    path batchkal

    output:
    path("${params.samplename}/matrix.abundance.mtx")
    
    storeDir params.outdir

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
  
    input:
    path index
    path reads

    output:
    path "${params.samplename}/output.bus"
    path "${params.samplename}/matrix.ec"
    path "${params.samplename}/transcripts.txt"
    
    storeDir params.outdir

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

/*
 * Step 3. Correct barcodes and sort
 * A future version of kallisto may include options to correct the UMI
 */
process corrsort {

    input:
    path outbus
    path white

    output:
    path("${params.samplename}/output.cor.sort.bus")
    
    storeDir params.outdir

    when: params.protocol!='plate'

    script:
    """
    mkdir -p ${params.samplename}
    bustools correct -w $white -o ${params.samplename}/output.cor.bus ${outbus}
    bustools sort -o ${params.samplename}/output.cor.sort.bus -t ${params.cores} ${params.samplename}/output.cor.bus
    """
}

/*
 * Step 4a. Obtain the UMI read counts
 */
process umicounts {

    input:
    path outbus

    output:
    path "${params.samplename}/umicount.txt"
    
    storeDir params.outdir

    when: params.protocol!='plate'

    script:
    """
    mkdir -p ${params.samplename}
    bustools text -o ${params.samplename}/umicount.txt ${outbus}
    """
}

/*
 * Step 4b. Obtain the counts
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

    when: params.protocol!='plate'

    script:
    """
    mkdir -p ${params.samplename}
    bustools count --em -t ${outtrans} -e ${outmat} -g $t2g --genecounts \\
    -o ${params.samplename}/genecounts ${outbus}
    """
}

/*
 * Step 5. Make Seurat object for plate
 */
process makeSeuratPlate {

    input:
    path outps
    path t2g

    output:
    path "${params.samplename}_srat.RDS"
    
    storeDir "${params.outdir}/${params.samplename}"

    when: params.protocol=='plate'

    script:
    """
    #!/usr/local/bin/Rscript --vanilla

    library(Seurat)

    # read in data
    topdir = "${params.samplename}" # source dir
    exp = Matrix::readMM("${outps}") #read matrix
    bc = read.csv(paste0(topdir, "/matrix.cells"), header = F, stringsAsFactors = F)
    g = read.csv(paste0(topdir, "/transcripts.txt"), header = F, stringsAsFactors = F)
    dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
    count.data = Matrix::t(exp)

    # summarise transcripts by gene name
    t2g = read.table("$t2g", header = F, stringsAsFactors = F, sep = "\t")
    exp_gene = rowsum(as.matrix(count.data), t2g\$V2)

    # create Seurat object
    srat = CreateSeuratObject(counts = exp_gene)

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
 * Step 5. Make Seurat object for 10x
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

    # read in data
    topdir = "${params.samplename}" # source dir
    exp = Matrix::readMM("${outmtx}") #read matrix
    bc = read.csv("${outbc}", header = F, stringsAsFactors = F)
    g = read.csv("${outg}", header = F, stringsAsFactors = F)
    dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
    count.data = Matrix::t(exp)

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
    umi = read.table("${umic}", sep = "\t", header = F, stringsAsFactors = F)
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
 * Step 5. Make Seurat object for ParseBio
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
    #!/usr/local/bin/Rscript --vanilla

    library(Seurat)
    library(DropletUtils)
    library(DelayedArray)

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
    #umi = read.table("${umic}", sep = "\t", header = F, stringsAsFactors = F)
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
 * Step Tiss. Get tissue alignment and fiducials
 */
process getTissue {
    storeDir "${params.outdir}/${params.samplename}"

    //input:
    //file imageal from file(params.imageal)
    //file imagef from file(params.imagef)

    output:
    path "mock/outs/spatial/aligned_fiducials.jpg"
    path "mock/outs/spatial/scalefactors_json.json"
    path "mock/outs/spatial/tissue_lowres_image.png"
    path "mock/outs/spatial/detected_tissue_image.jpg"
    path "mock/outs/spatial/tissue_hires_image.png"
    path "mock/outs/spatial/tissue_positions_list.csv"

    when: params.protocol=='visiumv1'

    script:
    """
    spaceranger count --id=mock --fastqs=./mock_fastq \\
    --transcriptome=/links/groups/treutlein/USERS/tomasgomes/gene_refs/human/refdata-gex-GRCh38-2020-A --image=${params.imagef} --slide=${params.images} --area=${params.imagear} --loupe-alignment=${params.imageal} --localcores=4
    """

}


/*
 * Step 5. Make Seurat object for Visium v1
 */
process makeSeuratVisium {

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

    when: params.protocol=='visiumv1'

    script:
    """
    #!/usr/local/bin/Rscript --vanilla

    library(Seurat)
    library(DropletUtils)

    # read in data
    topdir = "${params.samplename}" # source dir
    exp = Matrix::readMM("${outmtx}") #read matrix
    bc = read.csv(paste0(topdir, "/genecounts.barcodes.txt"), header = F, stringsAsFactors = F)
    g = read.csv(paste0(topdir, "/genecounts.genes.txt"), header = F, stringsAsFactors = F)
    dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
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

    # UMI duplication
    umi = read.table("${umic}", sep = "\t", header = F, stringsAsFactors = F)
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

    # create Seurat object
    srat <- CreateSeuratObject(counts = count.data, assay = "Spatial") # create object
    image <- Read10X_Image(image.dir = paste0(topdir, "/mock/outs/spatial/"),
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
 * Execute workflow
 */
workflow {
  starting()
  index(params.transcriptome)
  pseudoalPlate(index.out, batch_kal.collect())
  pseudoal(index.out, read_files_kallisto.collect())
  corrsort(pseudoal.out[0], bc_wl_kal.collect())
  umicounts(corrsort.out)
  countbus(corrsort.out, pseudoal.out[1], pseudoal.out[2], t2g_kal.collect())
  makeSeuratPlate(pseudoalPlate.out, t2g_plate.collect())
  makeSeurat10x(countbus.out[0], countbus.out[1], countbus.out[2], umicounts.out)
  makeSeuratParse(countbus.out[0], countbus.out[1], countbus.out[2], umicounts.out)
  //getTissue()
  //makeSeuratVisium(countbus.out[0], countbus.out[1], countbus.out[2], umicounts.out)
}
