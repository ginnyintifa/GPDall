

#' mapping vcf to various units
#'
#' @param vcf_folderPath the path to the folder of vcf files
#' @param grab_start_string a string indicating the sample's name
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
                        grab_start_string = NULL,
                        grab_sep = NULL,
                        grab_number = NULL,
                        mapping_vcf_to,
                        mapTo_fileName,
                        gtf_df,
                        reg_fileName,
                        ud_fileName,
                        output_folderPath,
                        output_tag
)

{
  # vcf_folderPath = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/vcf_folder_20210621/"
  # mapping_vcf_to = "protUnits"
  # mapTo_fileName = "/Users/ginny/Desktop/mapped_protUnits.tsv"
  # gtf_df = parse_gtf
  # reg_fileName = NULL
  # ud_fileName = NULL
  # output_folderPath = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/test_vcf_20210906/"
  # output_tag = "test_protUnits_sampleMapped"




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

      pb = patientInfo_extract(filenames[x],
                               grab_start_string = grab_start_string,
                               grab_sep = grab_sep,
                               grab_number = grab_number)

      matched = mapVCFtoProtUnits (vcf_file = filenames[x],
                                   protMappedGeno_file = mapTo_fileName,
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
        pb = patientInfo_extract(filenames[x],
                                 grab_start_string = grab_start_string,
                                 grab_sep = grab_sep,
                                 grab_number = grab_number)

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
          pb = patientInfo_extract(filenames[x],
                                   grab_start_string = grab_start_string,
                                   grab_sep = grab_sep,
                                   grab_number = grab_number)

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
            pb = patientInfo_extract(filenames[x],
                                     grab_start_string = grab_start_string,
                                     grab_sep = grab_sep,
                                     grab_number = grab_number)

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

  write.table(map_result, paste0(output_folderPath, output_tag,"_", "allMappedCount.tsv"),
              quote = F, row.names = F, sep = "\t")

if(nrow(map_result)>0)
{



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


  write.table(mat_df, paste0(output_folderPath, output_tag,"_", "countMatrix.tsv"),
              quote = F, row.names = F, sep = "\t")



  return(mat_df)
}else{
  cat("No variants mapped.","\n")
}

}





#
# gb = get_geneborder(gtf_df = parse_gtf,
#                 geneList = "TNFRSF4", # a list of gene symbols of interest
#                 geneBorder_filename = "test_reg.tsv")

