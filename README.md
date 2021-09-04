# GPDall

'gpd_workflow' is the main function in mapping variants from vcf file to various kinds of units. It will generate a output file with the detials of mapping in the 'output_folderPath',
and its return value is a matrix where each row is a unit and the number of mappable variants.

Please specify the type of units by the parameter 'mapping_vcf_to'. The following unit types are acceptable, "GTF", "regulatory", "protUnits" and "userDefine". With each type of unit specified, users are expected to specifiy the corresponding unit file in one of the following parameters, "mapTo_fileName" for "protUnits", "gtf_fileName" for "GTF", "reg_fileName" for "regulatory", and "ud_fileName" for "userDefine". 


For example, the following function maps vcf files in a fold to units seen in the parsed GTF file. 


```{r}
mat = gpd_workflow(vcf_folderPath = "/Path/to/vcf/file/folder/",
              mapping_vcf_to = "GTF",
              mapTo_fileName = NULL,
              gtf_fileName = "/path/to/parsed/gtf/file",
              reg_fileName = NULL, 
              ud_fileName = NULL,
              output_folderPath = "/path/to/your/output/folder",
              output_tag = "test_gtf")
```


