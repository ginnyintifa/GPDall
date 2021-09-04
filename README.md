# GPDall

## map to GTF regions from a GTF file 

'gpd_workflow' is the main function in mapping variants from vcf file to various kinds of units. It will generate a output file with the detials of mapping in the 'output_folderPath',
and its return value is a matrix where each row is a unit and the number of mappable variants.

Please specify the type of units by the parameter 'mapping_vcf_to'. The following unit types are acceptable "GTF", "regulatory", "protUnits" and "userDefine".


```{r}
mat = gpd_workflow(vcf_folderPath = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/vcf_folder_20210621/",
              mapping_vcf_to = "GTF",
              mapTo_fileName = NULL,
              gtf_fileName = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/parse_gtf.tsv",
              reg_fileName = NULL, 
              ud_fileName = NULL,
              output_folderPath = "/Users/ginny/Google Drive/R_GPD/GPD_package_0401/modiInput_202104/test_vcf_20210724/",
              output_tag = "test0728_gtf")
```
