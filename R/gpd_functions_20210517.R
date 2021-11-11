


#' Generate units based on UTR regions according to given downstream and upstream base pair numbers
#' @param up5TUR_bp number of base pairs upstream of 5' UTR, default to 1000
#' @param down3UTR_bp number of base pairs downstream of 3' UTR, default to 1000
#' @param gtf_border dataframe of genes' border info
#' @param regUnit_filename output filename of userdefined regions
#' @import dplyr data.table magrittr
#' @keywords get user defined regulatory regions from GTF file
#' @export

defineRegion_UTR = function(up5UTR_bp = 1000, ### how many base pair upstream of 5' UTR, default to 1000
                            down3UTR_bp = 1000, #### how many base pair downstream of 3' UTR, default to 1000
                            gtf_border,# = gene_boarder,
                            regUnit_filename)

{

  # #
  #   up5UTR_bp = 1000
  #    down3UTR_bp = 1000
  #     gtf_border = gb

  all_chros = unique(gtf_border$CHRO)

  get_unit = rbindlist(lapply(1:length(all_chros), function(x) {


   all_df = gtf_border[1,]%>%
     dplyr::mutate(up_start = NA, down_end = NA)

    this_chro = gtf_border%>%
      dplyr::filter(CHRO == all_chros[x], STRAND == "+")%>%
      dplyr::arrange(TSS)


    if(nrow(this_chro)>0)
    {
      up_start = unlist(lapply(1:nrow(this_chro), function(k) {

        it = max(1, this_chro$TSS[k]-up5UTR_bp)
        if(k>1)
        {
          itt = max(it, (this_chro$TES[k-1]+1))
        }else{
          itt = it
        }

        return(itt)
      }
      ))


      ####the downstream utr border of a gene should be smaller than the tss of the gene after it

      down_end = unlist(lapply(1:nrow(this_chro), function(k) {

        it =  this_chro$TES[k] + down3UTR_bp

        if(k<nrow(this_chro))
        {

          itt = min(it, (this_chro$TSS[k+1]-1))

        }else{
          itt = it
        }

        return(itt)

      }))


      this_df = this_chro%>%
        dplyr::mutate(up_start, down_end)


      all_df = this_df%>%
        na.omit()
    }



####the upper utr border of a gene should be biger than the tes of the gene in front of it


    #######  add another arm of the reverse strand thing
    ######  need to modify

    that_chro = gtf_border%>%
      dplyr::filter(CHRO == all_chros[x], STRAND == "-")%>%
      dplyr::arrange(desc(TSS))

 if(nrow(that_chro)>0)
 {
   rup_start = unlist(lapply(1:nrow(that_chro), function(k) {

     it =  that_chro$TSS[k] + up5UTR_bp
     if(k>1)
     {
       itt = min(it, (that_chro$TES[k-1]-1))
     }else{
       itt = it
     }

     return(itt)
   }
   ))



   rdown_end = unlist(lapply(1:nrow(that_chro), function(k) {

     it =  max(1,that_chro$TES[k] - down3UTR_bp)

     if(k<nrow(that_chro))
     {

       itt = max(it, (that_chro$TSS[k+1]+1))

     }else{
       itt = it
     }

     return(itt)

   }))


   that_df = that_chro%>%
     dplyr::mutate(up_start = rup_start,
                   down_end = rdown_end)

   all_df = rbind(all_df, that_df)%>%
     na.omit()


 }


    return(all_df)


  }))


  write.table(get_unit, regUnit_filename,
              quote = F, row.names = F, sep = "\t")


}


#' Generate prot units with genome coordinates
#' @param protUnit_filename file of protein units of interest
#' @param mappedProtUnit_filename output filename of protein units with genome coordinates added
#' @import dplyr data.table magrittr
#' @keywords get the border coordinates of a list of genes
#' @export
get_protGeno = function(protUnit_filename,
                        mappedProtUnit_filename)
{


  pe = fread(protUnit_filename,
             stringsAsFactors = F, data.table = F)

  p_g = rbindlist(lapply(1:nrow(pe), function(x) {

    #  if(x%%100 ==0)

    #  cat(x,"\n")
    po = IRanges(start = pe$start_position[x], end = pe$end_position[x], names = pe$uniprot_accession[x])

    r = proteinToGenome(po, EnsDb.Hsapiens.v86, idType = "uniprot_id")

    ##### gather all of them and take the union of chromosome and position

    if(length(r[[1]])>0)
    {
      if("unlistData" %in% slotNames(r[[1]]))
      {
        df = r[[1]]@unlistData

        CHROM = rep(as.character(df@seqnames@values[1]), length(df))
        gstart = df@ranges@start
        strand = as.character(df@strand@values[1])

        gend = df@ranges@width + df@ranges@start -1


        this_df = data.frame(CHROM, gstart, gend,strand,stringsAsFactors = F)%>%
          unique()


      }else{
        df = r[[1]]

        CHROM = rep(as.character(df@seqnames@values[1]), length(df))
        gstart = df@ranges@start
        strand = as.character(df@strand@values[1])
        gend = df@ranges@width + df@ranges@start -1


        this_df = data.frame(CHROM, gstart, gend, strand, stringsAsFactors = F)%>%
          unique()

      }

      lb = unique(this_df)


      final_df = cbind(pe[rep(x,nrow(lb)),], lb)

      return(final_df)

    }


  }))


  write.table(p_g, mappedProtUnit_filename,
              quote = F, row.names = F, sep = "\t")

}





#' Generate gene borders according to the gtf file given a list of genes of interest
#' @param gtf_df parsed gtf datafrome, can be loaded from the package
#' @param geneList a list of gene symbols of interest
#' @param geneBorder_filename output filename of borders of genes of interest
#' @import dplyr data.table magrittr
#' @keywords get ghe border coordinates of a list of genes
#' @export

get_geneborder = function(gtf_df,
                          geneList,
                          geneBorder_filename)
{

  # gtf_df = parse_gtf
  # geneList = gl
  # geneBorder_filename = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/test_vcf_20210906/geneBorders.tsv"


  real_genes = gtf_df%>%
    dplyr::filter(transcript_type == "protein_coding", gene_name %in% geneList)

  gene_sets = unique(real_genes$gene_name)
  gene_border = rbindlist(lapply(1:length(gene_sets), function(x) {
    get_gene = gtf_df%>%
      dplyr::filter(gene_name == gene_sets[x])

    tss = min(get_gene$START)
    tes = max(get_gene$END)

    if(get_gene$STRAND[1] == "-"){
      tss = max(get_gene$END)
      tes = min(get_gene$START)

    }


    df = data.frame(gene_name = gene_sets[x],gene_id = get_gene$gene_id[1],
                    CHRO = get_gene$CHRO[1], STRAND = get_gene$STRAND[1],
                    TSS = tss, TES = tes, stringsAsFactors = F)


    return(df)


  }))


  gene_border = gene_border%>%
    dplyr::arrange(CHRO, TSS)

  write.table(gene_border, geneBorder_filename,
              quote = F, row.names = F, sep = "\t")

  return(gene_border)

}






mapVCFtoProtUnits = function(vcf_file,
                             protMappedGeno_file,
                             vcfMapped_prot_outputName,
                             vcfMapped_protUnitCount_outputName)
{

  vcf = fread(vcf_file, skip = "#CHROM",
              stringsAsFactors = F, data.table = F)

  colnames(vcf) = gsub("#","",colnames(vcf))


 # if(genoMapped == T)
  #{
    p_g = fread(protMappedGeno_file, stringsAsFactors = F, data.table = F)

  #}else{

  #   pe = fread(protUnit_file,
  #              stringsAsFactors = F, data.table = F)
  #
  #   p_g = rbindlist(lapply(1:nrow(pe), function(x) {
  #
  #     if(x%%100 ==0)
  #       cat(x,"\n")
  #     po = IRanges(start = pe$start_position[x], end = pe$end_position[x], names = pe$uniprot_accession[x])
  #
  #     r = proteinToGenome(po, EnsDb.Hsapiens.v86, idType = "uniprot_id")
  #
  #     ##### gather all of them and take the union of chromosome and position
  #
  #     if(length(r[[1]])>0)
  #     {
  #       if("unlistData" %in% slotNames(r[[1]]))
  #       {
  #         df = r[[1]]@unlistData
  #
  #         CHROM = rep(as.character(df@seqnames@values[1]), length(df))
  #         gstart = df@ranges@start
  #         strand = as.character(df@strand@values[1])
  #
  #         gend = df@ranges@width + df@ranges@start -1
  #
  #
  #         this_df = data.frame(CHROM, gstart, gend,strand,stringsAsFactors = F)%>%
  #           unique()
  #
  #
  #       }else{
  #         df = r[[1]]
  #
  #         CHROM = rep(as.character(df@seqnames@values[1]), length(df))
  #         gstart = df@ranges@start
  #         strand = as.character(df@strand@values[1])
  #         gend = df@ranges@width + df@ranges@start -1
  #
  #
  #         this_df = data.frame(CHROM, gstart, gend, strand, stringsAsFactors = F)%>%
  #           unique()
  #
  #       }
  #
  #       lb = unique(this_df)
  #
  #
  #       final_df = cbind(pe[rep(x,nrow(lb)),], lb)
  #
  #       return(final_df)
  #
  #     }
  #
  #
  #   }))
  #
  #
  #   write.table(p_g, protMappedGeno_outputName,
  #               quote = F, row.names = F, sep = "\t")
  #
  # }

  matched = rbindlist(lapply(1:nrow(vcf), function(x) {


    #cat(x, "\n")
    this_chrom = gsub("chr","", vcf$CHROM[x])
    this_pos = vcf$POS[x]


    get_p_g =  p_g%>%
      dplyr::filter(CHROM == this_chrom)%>%
      dplyr::filter(gstart <= this_pos, gend >= this_pos)

    if(nrow(get_p_g)>0)
    {

      p_g_info = get_p_g%>%
        dplyr::mutate(prot_info = paste(uniprot_accession, gene_name, start_position, end_position, unit_name, unit_label, sep = "_"),
                      geno_info = paste(CHROM, gstart, gend, strand, sep = "_"))%>%
        dplyr::select(prot_info, geno_info)

      this_df = cbind(p_g_info, vcf[rep(x, nrow(p_g_info)),])


      return(this_df)
    }



  }),use.names=TRUE)


  write.table(matched, vcfMapped_prot_outputName,
              quote = F, row.names = F, sep = "\t")

  cat("protein units mapped to genome. ","\n")


  if(nrow(matched)>0)
  {

    matched_count = matched%>%
      dplyr::select(-geno_info)%>%
      dplyr::group_by(prot_info)%>%
      dplyr::mutate(count = n())



    write.table(matched_count, vcfMapped_protUnitCount_outputName,
                quote = F, row.names = F, sep = "\t")

   return(matched_count)

  }
}


mapVCFtoGTF= function(vcf_file,
                      gtf,   ###change to parsed gtf file 20210709
                      vcfMapped_gtf_outputName,
                      vcfMapped_gtfCount_outputName,
                      vcfMapped_gtf_nonTranslated_outputName,
                      vcfMapped_gtf_nonTranslatedCount_outputName)

{



  match_gtf = function(x,vcf,gtf)
  {

    this_chro = vcf$CHROM[x]
    this_pos = vcf$POS[x]

    fil_gtf =gtf%>%
      dplyr::filter(CHRO == this_chro &START <= this_pos & END >= this_pos)

    if(nrow(fil_gtf)>0)
    {


      parse_gtf = fil_gtf%>%
        unique()%>%
        dplyr::mutate(vcf_chro = this_chro,
                      vcf_pos = this_pos)%>%
        dplyr::select(vcf_chro, vcf_pos, everything())


    }else{

      parse_gtf = data.frame(vcf_chro = this_chro,
                             vcf_pos = this_pos,
                             CHRO = this_chro,
                             SOURCE = "",
                             FEATURE = "not_specified",
                             START = "",
                             END = "",
                             STRAND = "",
                             PHASE = "",
                             gene_id = "",gene_type = "", gene_name = "",
                             transcript_id  = "", transcript_type = "", transcript_name = "",
                             exon_number = "", exon_id ="")


    }

    return(parse_gtf)

  }




  vcf = fread(vcf_file, skip = "#CHROM",
              stringsAsFactors = F, data.table = F)
  colnames(vcf) = gsub("#","",colnames(vcf))

  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)


  vcf_mapped_gtf = rbindlist(foreach(x= 1:nrow(vcf), .packages = "dplyr")%dopar% match_gtf(x,vcf,gtf))



  stopCluster(cl)


  ################################
  ###############################
  ###############################
  ##############################



  vcf_mapped_gtf = vcf_mapped_gtf%>%
    dplyr::mutate(gtf_info = paste(CHRO, SOURCE, FEATURE, START, END, STRAND, PHASE, gene_name, transcript_name, sep = "_"))%>%
    unique()

  write.table(vcf_mapped_gtf, vcfMapped_gtf_outputName,
              quote = F, row.names = F, sep  = "\t")


  vcf_mapped_gtf_counts = vcf_mapped_gtf%>%
    dplyr::select(gtf_info,CHRO, SOURCE, FEATURE, START, END, STRAND, PHASE, gene_name, transcript_name)%>%
    dplyr::group_by(gtf_info)%>%
    dplyr::mutate(count = n())%>%
    unique()



  write.table(vcf_mapped_gtf_counts, vcfMapped_gtfCount_outputName,
              quote = F, row.names = F, sep  = "\t")



  vcf_mapped_gtf_nonTranslated = vcf_mapped_gtf%>%
    dplyr::filter(FEATURE == "UTR" | FEATURE == "not_specified")%>%
    unique()


  write.table(vcf_mapped_gtf_nonTranslated, vcfMapped_gtf_nonTranslated_outputName,
              quote = F, row.names = F, sep  = "\t")


  vcf_mapped_gtf_nonTranslated_counts = vcf_mapped_gtf_nonTranslated%>%
    dplyr::select(gtf_info,CHRO, SOURCE, FEATURE, START, END, STRAND, PHASE, gene_name, transcript_name)%>%
    dplyr::group_by(gtf_info)%>%
    dplyr::mutate(count = n())%>%
    unique()

  write.table(vcf_mapped_gtf_nonTranslated_counts, vcfMapped_gtf_nonTranslatedCount_outputName,
              quote = F, row.names = F, sep  = "\t")
  return(vcf_mapped_gtf_counts)

}




#### this can be parallezed too

mapVCFtoReg = function(vcf_file,
                       reg_file,
                       vcfMapped_reg_outputName,
                       vcfMapped_regCount_outputName)


{


  vcf = fread(vcf_file, skip = "#CHROM",
              stringsAsFactors = F, data.table = F)
  colnames(vcf) = gsub("#","",colnames(vcf))

  reg =  fread(reg_file,
               stringsAsFactors = F, data.table = F)



  match_reg = function(x,vcf,reg)
  {

    this_chro = vcf$CHROM[x]
    this_pos = vcf$POS[x]

    fil_reg_up =reg%>%
      dplyr::filter(CHRO == this_chro & (up_start <= this_pos & TSS >= this_pos))%>%
      dplyr::mutate(reg_spec = "upstreamTSS")%>%
      dplyr::mutate(vcf_chro = this_chro,
                    vcf_pos = this_pos)

    fil_reg_down =reg%>%
      dplyr::filter(CHRO == this_chro & (TES <= this_pos & down_end >= this_pos))%>%
      dplyr::mutate(reg_spec = "downstreamTES")%>%
      dplyr::mutate(vcf_chro = this_chro,
                    vcf_pos = this_pos)

    fil_reg = rbind(fil_reg_up, fil_reg_down)%>%
      dplyr::select(vcf_chro, vcf_pos, everything())

    ##### parse this portion only

    if(nrow(fil_reg)>0)
    {
      fil_reg = fil_reg
    }else{


      fil_reg = data.frame(vcf_chro = this_chro,
                           vcf_pos = this_pos,
                           gene_name = "",
                           gene_id = "",
                           CHRO = "",
                           STRAND = "",
                           TSS = "",
                           TES = "",
                           up_start = "",
                           down_end = "",
                           reg_spec = "",
                           stringsAsFactors = F)


    }



    return(fil_reg)

  }

  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)


  vcf_mapped_reg = rbindlist(foreach(x= 1:nrow(vcf), .packages = "dplyr")%dopar% match_reg(x,vcf,reg))



  stopCluster(cl)


  ###############################


    vcf_mapped_reg = vcf_mapped_reg%>%
      dplyr::mutate(reg_info = paste(CHRO, STRAND,gene_name, up_start, TSS, TES, down_end,reg_spec, sep = "_"))%>%
      unique()%>%
    dplyr::filter(nchar(reg_info)>7)

    write.table(vcf_mapped_reg, vcfMapped_reg_outputName,
                quote = F, row.names = F, sep  = "\t")

    vcf_mapped_reg_counts = vcf_mapped_reg%>%
      dplyr::select(reg_info,CHRO, STRAND, gene_name, gene_id, TSS, TES, up_start, down_end, reg_spec)%>%
      dplyr::group_by(reg_info)%>%
      dplyr::mutate(count = n())%>%
      unique()

    write.table(vcf_mapped_reg_counts, vcfMapped_regCount_outputName,
                quote = F, row.names = F, sep  = "\t")

    return(vcf_mapped_reg_counts)




}





mapVCFtoUD = function(vcf_file,
                      ud_file,
                      vcfMapped_ud_outputName,
                      vcfMapped_udCount_outputName)


{

  vcf = fread(vcf_file, skip = "#CHROM",
              stringsAsFactors = F, data.table = F)
  colnames(vcf) = gsub("#","",colnames(vcf))

  ud =  fread(ud_file,
              stringsAsFactors = F, data.table = F)

  d_ud = rbindlist(lapply(1:nrow(vcf), function(x) {

    this_chro = vcf$CHROM[x]
    this_pos = vcf$POS[x]

    fil_ud =ud%>%
      dplyr::filter(CHRO == this_chro & (start <= this_pos & end >= this_pos))%>%
      dplyr::mutate(vcf_chro = this_chro,
                    vcf_pos = this_pos)


    if(nrow(fil_ud)>0)
    {


      return(fil_ud)


    }



  }))


  if(nrow(d_ud)>0)
  {

    vcf_mapped_ud = d_ud%>%
      dplyr::mutate(ud_info = paste(CHRO, STRAND,gene_name,start, end, sep = "_"))%>%
      unique()

    write.table(vcf_mapped_ud, vcfMapped_ud_outputName,
                quote = F, row.names = F, sep  = "\t")



    vcf_mapped_ud_counts = vcf_mapped_ud%>%
      dplyr::select(ud_info,CHRO, STRAND, gene_name, gene_id, start, end)%>%
      dplyr::group_by(ud_info)%>%
      dplyr::mutate(count = n())%>%
      unique()



    write.table(vcf_mapped_ud_counts, vcfMapped_udCount_outputName,
                quote = F, row.names = F, sep  = "\t")

  }

  return(vcf_mapped_ud_counts)
}




patientInfo_extract =function(vcf_file,
                              grab_start_string,
                              grab_sep,
                              grab_number){



  f = readLines(vcf_file)
  f_tf = grepl("##TUMOR", f)
  if(sum(f_tf)==1)
  {
    f_t = grep("##TUMOR", f, value = T)
    f_t_split = unlist(strsplit(f_t, split = ","))
    b = unlist(strsplit(f_t_split[1],split = "Sample="))[2]
  }else{

    seps = unlist(strsplit(vcf_file, split = grab_sep, fixed = T))

    w = which(seps == grab_start_string)

    wend = w+grab_number-1

    tcga_tag = paste(seps[w:wend], collapse = grab_sep)

    b = tcga_tag
  }



  return(b)
}



