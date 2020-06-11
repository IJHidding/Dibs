# Dibs
The DIBS function performs imputation benchmarking for datasets containing NA's or zeros. 

## Installing
To install and use this library in your R, download DIBS.1.0.0.tar.gz and then run:

    install.packages(path_to_file, repos = NULL, type="source")

Where path_to_file will be something like "C:\\DIBS_1.0.0.tar.gz" on windows.
Or "/home/DIBS_1.0.0.tar.gz" on UNIX. 

Package: DIBS
Type: Package
Title: Data Imputation Benchmark Software
Version: 1.0.0
Author: Iwan Hidding
Maintainer: iwan_hidding@live.nl
Description: The DIBS function performs imputation benchmarking for datasets containing NA's or zeros. 
    
    Usage
    dibs <- (x, ...)
    
    ## Default methods
    dibs <- function(x, transform_method="none",variable_on="column", impute_setting=1, impute_method, 
    listofzeros=c(10, 25, 50, 75, 90), imp_method=c("norm", "pmm", "rf", "knn"), cor_method="pearson", 
    numberofmetabolites= 1, cutoff = 0.55, minimum_metabolites=4, maximum_metabolites=10, low_cor_cutoff=0.4) 
    
    Arguments
    
    x                   A dataframe containing zeros or NA's to be used in the benchmark or actual imputation. 
    
    transform_method    The method to be used for transforming or scaling the data; options available are:"none", 
                        "scale", "minmax", "log", "zscore"
    
    variable_on         An indication of the structure of the input data. Indicate if the samples are on row or column,
                        options available: "row", "column"
    
    impute_setting      Impute for the function, available options: 1: Benchmark the function to return the best imputation 
                        method, setting 2: benchmarks the best method and then applies this method to impute the dataset, 
                        returns the imputed dataset. Setting 3: needs the impute_method variable and only performs the imputation
                        without benchmarking first. For when the best  method is already known.
    
    impute_method       The imputation method to use for imputation when impute_setting=3 is used. 
    
    listofzeros         The different percentages of missingness to be tested for the benchmarking. A wide range is set as default
                        it is recommended to adapt this to your own dataset.
    
    imp_method          The imputation methods available in this function, options are the Knn method and all methods in mice.impute,
                        for an overview use methods(mice) after loading the mice library. The "rf", "norm" and "pmm" methods have been 
                        tested and are supported.
    
    cor_method          The correlation method used to determine the correlation between the variables. With "pearson" as default
                        and "kendall", or "spearman" as other options. 
    
    numberofmetabolites The percentage of variables to be tested in the benchmarking setup. Testing has shown around 75 variables are
                        needed to get a consistent recommendation. Default is 1%, adapt this to reach at least 75 variables in your own
                        dataset. 
    
    cutoff              The minimal required correlation score for other variables to count as correlated for imputation.
    
    minimum_metabolites The minimal number of other variables to be used for the benchmarking.
    
    maximum_metabolites The maximum number of other variables to be used for the benchmarking if there are many correlated variables. 
                        Recommended to keep below 25.
    
    low_cur_cutoff      The minimum correlation score for the imputation to be considered succesful. Any variables scoring below here are
                        recorded and stored in the output. 
                        
    
    Details
    
    It will introduce NA's in data where it does not yet exist, then impute those values back in using variables 
    with a high correlation to the variable to be imputed. It compares the new values with the original values,
    then it averages all those results to detect the overal best method. Depending on the mode selected, DIBS will 
    then automatically select this method and impute your imput dataset. With setting 2 or 3 it will return the input data
    frame with filled in zeros. If setting 1 is selected DIBS returns the correlation dataframe. Additionally an output pdf
    and csv file will be created. The csv contains the correlation dataframe for review and the pdf contains several plots;
    a histogram of the overall missingness; a heatmap showing correlation between random variables from the dataset, giving an
    indication if there is any correlation between the variables; density plots showing a percentage of well performing and 
    poorly performing variables in the dataset. A list of poorly performing variables used in the benchmark is added at the 
    end aswell as any variables for which not enough correlated other variables could be found. 
Depends: R (>= 3.6.0)

Imports:
    dplyr,
    ggplot2,
    grid,
    gridExtra,
    mice,
    reshape2
License: GNU Lesser General Public License
Encoding: UTF-8
LazyData: true

# Next steps
This is version 1.0.0 of the DIBS tool. It has been created based on metabolite data, but its been made available to be used on other types
of omics data. As stated above the tested methods are the norm, pmm and rf methods from the mice package. Other methods can be used but have
to be tested to ensure they work properly. Additionally other types of data should be tested with this setup to ensure proper performance and 
behaviour. 

Performance is an additional big part. Multiprocessing options should be relatively easy to implement and will significantly enhance performance,
as the whole benchmark can be split up by, for example, metabolite and this can be split over all available cores which will be likely many times 
faster than the current runtime. 
