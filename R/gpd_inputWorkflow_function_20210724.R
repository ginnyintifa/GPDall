

#' mapping vcf to various units
#'
#' @param vcf_folderPath the path to the folder of vcf files
#' @param mapping_vcf_to mapping mode, one of "GTF","regulatory","protUnits","userDefine"
#' @param mapTo_fileName file of protein units to be mapped with genome coordinates
#' @param gtf_df gtf dataframe to be mapped, can be loaded with the package
#' @param reg_fileName file of regulatory units to be mapped
#' @param ud_fileName file of user-defined regions
#' @param output_folderPaht the directory holding output files
#' @param output_tag a tag for output files
#' @import dplyr data.table magrittr doParallel foreach
#' @keywords mapping vcf to various units
#' @export
#' @examples

gpd_workflow = function(vcf_folderPath,
                        mapping_vcf_to,
                        mapTo_fileName,
                        gtf_df,
                        reg_fileName,
                        ud_fileName,
                        output_folderPath,
                        output_tag
)

{


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
                                  # genoMapped = T,  ### new version change it into two steps
                                   protMappedGeno_file = mapTo_fileName,
                                  # protUnit_file= NULL,
                                   #protMappedGeno_outputName = NULL,
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
      gtf_df = gtf_df

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
