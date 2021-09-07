# GPDall

'gpd_workflow' is the main function in mapping variants from vcf file to various kinds of units. It will generate a output file with the detials of mapping in the 'output_folderPath', and its return value is a matrix where each row is a unit and the number of mappable variants.

'vcf_folderPath' is the path to the folder that holds the vcf files, each file is for an individual. 

Please specify the type of units by the parameter 'mapping_vcf_to'. The following unit types are acceptable, "GTF", "regulatory", "protUnits" and "userDefine". With each type of unit specified, users are expected to specifiy the corresponding unit file in one of the following parameters, "mapTo_fileName" for "protUnits", "gtf_fileName" for "GTF", "reg_fileName" for "regulatory", and "ud_fileName" for "userDefine". 


An annotated version of gtf file is included in the package, users can call it with `parse_gtf`, it is generated from Gencode "gencode.v38.annotation.gtf". https://www.gencodegenes.org/human/.


## Map variants to GTF file 


For example, the following function maps vcf files in a folder to units seen in the parsed GTF file. 

```{r}
mat = gpd_workflow(vcf_folderPath = "/Path/to/vcf/file/folder/",
              mapping_vcf_to = "GTF",
              mapTo_fileName = NULL,
              gtf_df = parse_gtf,
              reg_fileName = NULL, 
              ud_fileName = NULL,
              output_folderPath = "/path/to/your/output/folder",
              output_tag = "test_gtf")
```

## Map variants to protein units 

User can view a sample of the input protUnit_file by calling "protUnit_example"


```{r}
mat = gpd_workflow(vcf_folderPath = "/Path/to/vcf/file/folder/",
              mapping_vcf_to = "protUnits",
              mapTo_fileName = "/Path/to/protUnit_file",
              gtf_df = parse_gtf,
              reg_fileName = NULL, 
              ud_fileName = NULL,
              output_folderPath = "/path/to/your/output/folder",
              output_tag = "test_protUnits")
```


Under "protUnits" mode, if the genome coordinates are unknown, user can use the function `get_protGeno` to preprocess the protein unit information file. Please install two bioconductor packages as follows:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnsDb.Hsapiens.v86")

BiocManager::install("ensembldb")

```

And add genome coordinates to the protein unis in the following way:

```{r}

get_protGeno(protUnit_filename = unMapped_protUnit_example,
             mappedProtUnit_filename = "mapped_filename")

```


## Map variants to user-defined genome regions 

User can view a sample of the input userDefine_file by calling "userDefine_example"

```{r}
mat = gpd_workflow(vcf_folderPath = "/Path/to/vcf/file/folder/",
              mapping_vcf_to = "userDefine"
              mapTo_fileName = NULL,
              gtf_df = parse_gtf,
              reg_fileName = NULL, 
              ud_fileName = "/Path/to/userDefine_file",
              output_folderPath = "/path/to/your/output/folder",
              output_tag = "test_userDefine")
```


## Map variants to UTR based regulatory regions with user-defined downstream and upstream length 

User can view a sample of the input regulatory_file by calling "regulatoryRegion_example"

```{r}
mat = gpd_workflow(vcf_folderPath = "/Path/to/vcf/file/folder/",
              mapping_vcf_to = "regulatory"
              mapTo_fileName = NULL,
              gtf_df = parse_gtf,
              reg_fileName = NULL, 
              ud_fileName = "/Path/to/regulatory_file",
              output_folderPath = "/path/to/your/output/folder",
              output_tag = "test_regulatory")
```

Under "regulatory" mode, user may obtain the desired regulatory unit file by first obtaining the borders of genes of interest and then obtaining the downstream and upstream UTR regions by the follwing two functions 

```{r}


gtf_border = function(gtf_df,
	                  geneList, # a list of gene symbols of interest
                      geneBorder_filename)



defineRegion_UTR = function(up5UTR_bp = 1000, ### how many base pair upstream of 5' UTR, default to 1000
                            down3UTR_bp = 1000, #### how many base pair downstream of 3' UTR, default to 1000
                            gtf_border,
                            regUnit_filename)

```





