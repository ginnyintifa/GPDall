





#' Compare between 2 conditions (paired) in 1D
#'
#' @param vcf_folderPath names of the columns in the first condition
#' @param mapping_vcf_to names of the columns in the second condition
#' @param mapTo_fileName dataframe for the data, rows are genes, columns are patients in 2 conditions
#' @param gtf_fileName genes with values in at least the number of samples are included in comparison
#' @param reg_fileName seeds for the permutation
#' @param ud_fileName number of permutations
#' @param output_folderPaht the directory output files will be deposited in
#' @param output_tag a label for the anlaysis
#' @import dplyr data.table magrittr doParallel foreach
#' @keywords mapping
#' @export
#' @examples

gpd_workflow = function(vcf_folderPath,
                        mapping_vcf_to,
                        mapTo_fileName,
                        gtf_fileName,
                        reg_fileName,
                        ud_fileName,
                        output_folderPath,
                        output_tag
)

{
  #

  #vcf_folderPath = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/vcf_folder_20210621/"
  # mapping_vcf_to = "protUnits"
  # mapTo_fileName = "/Users/Ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/protGenoMapped.tsv"
  # gtf_fileName = NULL
  # reg_fileName = NULL
  # ud_fileName = NULL
  # output_folderPath = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/test_vcf_20210724/"
  # output_tag = "test0728_protUnits"
  #


  files = list.files(vcf_folderPath)
  filenames = paste(vcf_folderPath, files, sep = "")

  map_result = NULL

  ##################################################################
  ######################################################### MAP TO PROTEIN UNITS
  ##################################################################

  if(mapping_vcf_to == "protUnits")
  {

    map_result = rbindlist(lapply(1:length(filenames), function(x) {
      cat(x, "\n")

      pb = patientInfo_extract(filenames[x])

      matched = mapVCFtoProtUnits (vcf_file = filenames[x],
                                   genoMapped = T,
                                   protMappedGeno_file = mapTo_fileName,
                                   protUnit_file= NULL,
                                   protMappedGeno_outputName = NULL,
                                   vcfMapped_prot_outputName = paste0(output_folderPath, output_tag,"_",
                                                                      pb,"_", "mapToProt.tsv"),
                                   vcfMapped_protUnitCount_outputName = paste0(output_folderPath, output_tag,"_",
                                                                               pb,"_", "mapToProtCount.tsv"))

      if(is.null(matched)==F)
      {
        m_r = matched%>%
          dplyr::mutate(barcode = pb)

        return(m_r)
      }
    }), use.names = T)


  }else
  {
    ##################################################################
    ######################################################### MAP TO GTF file
    ##################################################################
    if(mapping_vcf_to == "GTF")
    {
      gtf_df = fread(gtf_fileName, stringsAsFactors = F, data.table = F)

      #### a series of operations
      map_result = rbindlist(lapply(1:length(filenames), function(x) {

        cat(x, "\n")
        pb = patientInfo_extract(filenames[x])

        tt = Sys.time()
        matched = mapVCFtoGTF(vcf_file = filenames[x],
                              gtf = gtf_df,
                              vcfMapped_gtf_outputName =  paste0(output_folderPath, output_tag,"_",
                                                                 pb,"_", "mapToGTF.tsv"),
                              vcfMapped_gtfCount_outputName =  paste0(output_folderPath, output_tag,"_",
                                                                      pb,"_", "mapToGTFCount.tsv"),
                              vcfMapped_gtf_nonTranslated_outputName =  paste0(output_folderPath, output_tag,"_",
                                                                               pb,"_", "mapToGTF_nonTranslated.tsv"),
                              vcfMapped_gtf_nonTranslatedCount_outputName =  paste0(output_folderPath, output_tag,"_",
                                                                                    pb,"_", "mapToGTFCount_nonTranslated.tsv"))

        et = Sys.time()

        if(is.null(matched)==F)
        {
          m_r = matched%>%
            dplyr::mutate(barcode = pb)

          return(m_r)
        }
      }), use.names = T)
    }else{
      ##################################################################
      ######################################################### MAP TO REGULATORY REGIONS
      ##################################################################

      if(mapping_vcf_to == "regulatory")
      {

        map_result = rbindlist(lapply(1:length(filenames), function(x) {
          cat(x, "\n")
          pb = patientInfo_extract(filenames[x])

          matched = mapVCFtoReg(vcf_file = filenames[x],
                                reg_file = reg_fileName,
                                vcfMapped_reg_outputName = paste0(output_folderPath, output_tag,"_",
                                                                  pb,"_", "mapToReg.tsv"),
                                vcfMapped_regCount_outputName = paste0(output_folderPath, output_tag,"_",
                                                                       pb,"_", "mapToRegCount.tsv"))

          if(is.null(matched)==F)
          {
            m_r = matched%>%
              dplyr::mutate(barcode = pb)

            return(m_r)
          }
        }), use.names = T)
      }else{
        ##################################################################
        ######################################################### MAP TO USER DEFINED REGIONS
        ##################################################################

        if(mapping_vcf_to == "userDefine")
        {
          map_result = rbindlist(lapply(1:length(filenames), function(x) {
            cat(x, "\n")
            pb = patientInfo_extract(filenames[x])

            matched = mapVCFtoUD(vcf_file = filenames[x],   #### can be parallelized too, advance this later
                                 ud_file = ud_fileName,
                                 vcfMapped_ud_outputName = paste0(output_folderPath, output_tag,"_",
                                                                  pb,"_", "mapToUD.tsv"),
                                 vcfMapped_udCount_outputName = paste0(output_folderPath, output_tag,"_",
                                                                       pb,"_", "mapToUDCount.tsv"))

            if(is.null(matched)==F)
            {
              m_r = matched%>%
                dplyr::mutate(barcode = pb)

              return(m_r)
            }

          }),use.names = T)

        }else{
          cat( "wrong mapping type")

        }
      }
    }
  }

  ######### before the final return need to make it a matrix

  write.table(map_result, paste0(output_folderPath, output_tag,"_", "allMappedCount.tsv"))


  #### I require what ever the type is the first column contains the unit info
  colnames(map_result)[1] = "unit_info"


  unique_units = unique(map_result$unit_info)
  patients = unlist(lapply(1:length(filenames), function(x) patientInfo_extract(filenames[x])))

  mat = matrix(0,length(unique_units), length(patients))

  for(i in 1:ncol(mat))
  {

    this_p = patients[i]

    sli = map_result%>%
      dplyr::filter(barcode == this_p)

    ur = match(sli$unit_info, unique_units)

    mat[ur,i] = sli$count

  }

  rs = rowSums(mat)

  mat_df = data.frame(unique_units, mat, rs,stringsAsFactors = F)
  colnames(mat_df) = c("unit_info", patients,"sum")

  mat_df = mat_df%>%
    dplyr::arrange(desc(sum))


  write.table(mat_df, paste0(output_folderPath, output_tag,"_", "countMatrix.tsv"))



  return(mat_df)

}
