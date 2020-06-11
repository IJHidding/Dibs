# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(mice)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(grid)
dibs <- function(dataset, transform_method="none",variable_on="column", impute_setting=1, impute_method, listofzeros=c(10, 25, 50, 75, 90), imp_method=c("norm", "pmm", "rf", "knn"), cor_method="pearson", numberofmetabolites= 1, cutoff = 0.55, minimum_metabolites=4, maximum_metabolites=10, low_cor_cutoff=0.4) {
  # A quick check to ensure that the data is in the right format for the rest of the analysis
  if (variable_on == "row") {
    dataset <- as.data.frame(t(dataset))
  }
  # Open the output file
  pdf(file="Benchmark_output.pdf",
      width = 10,
      height = 8)
  # Determine missingness and plot a histogram of the missingness
  missingvector <- apply(dataset, 1, function (x) { sum(x==0) / length(x) })
  hist(missingvector, main= "Missingness overview whole dataset")

  # Make the heatmap
  cormat <- round(cor(dataset[sample(nrow(dataset), 30), sample(ncol(dataset), 30)]),2)

  melted_cormat <- melt(cormat)
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    coord_fixed() + ggtitle("Heatmap of correlation of random subset of variables.")
  # Print the heatmap
  plot(ggheatmap)

  # Impute settings one and two perform the benchmark, this is checked for here.
  if (impute_setting == 1 || impute_setting == 2) {
    # Creates and empty dataframe of all the rows and columns in the output to be filled in
    df <- data.frame(matrix(ncol = length(imp_method), nrow=round(numberofmetabolites * 0.01 * ncol(dataset), 0) * length(listofzeros)))
    colnames(df) <- c(imp_method)

    listofmetabolites <- colnames(dataset)[c(sample(1: ncol(dataset), round((numberofmetabolites * 0.01 * ncol(dataset)), 0)))]
    zeroMetalist <- c()
    # This part adds the percentages of zeros tested to the metabolite name to keep track if performance changes
    for ( numberofzeros in listofzeros) {
      zeroMetalist <- c(zeroMetalist,  paste(listofmetabolites, paste(numberofzeros, "%", sep=""), sep="_"))
    }

    rownames(df) <- zeroMetalist

    missingness_count <- 0
    number_of_skipped_metabolites <- 0
    low_cor_vector <<- c()
    # first loops over all different zeros and then all different metabolites
    # This way all metabolites will have the same samples removed for each number of zeros
    for ( numberofzeros in listofzeros) {
      listofskippedmetabolites <- c()
      zeros <- sample(1:nrow(dataset), (numberofzeros * 0.01 * nrow(dataset)))
      for( metabolite in listofmetabolites) {
        zerometa <- paste(metabolite, paste("_", numberofzeros, "%", sep=""), sep="")
        print(paste(metabolite, paste("_", (numberofzeros * 0.01 * nrow(dataset)), sep=""), sep=""))
        new_x <- dataset
        saved_data <- new_x[c(zeros), metabolite]
        plotting_df <- as.data.frame(saved_data)
        new_x[c(zeros),metabolite] <- NA
        cor_df <- new_x[complete.cases(new_x),]
        cor_matrix <- as.data.frame(cor(cor_df, method = cor_method))
        sorteddata <- cor_matrix[order(-cor_matrix[metabolite]),]
        number_of_correlated_metabolites <- length(na.omit(rownames((sorteddata[which(sorteddata[,metabolite]> cutoff),]))[1:maximum_metabolites]))
        # A quick check if there are enough correlated metabolites available for a proper imputation, otherwise this variable is skipped and noted
        if (number_of_correlated_metabolites < minimum_metabolites ){
          listofskippedmetabolites <- c(listofskippedmetabolites, metabolite)
          number_of_skipped_metabolites <- number_of_skipped_metabolites + 1
          next
        }
        relevantstuff <- rownames(sorteddata)[1:number_of_correlated_metabolites]
        mice_dataframe <- (select(new_x, relevantstuff))

        # The only method that is not available in mice that is available in dibs is the knn,
        # so if called it will be skipped here and called in the next step
        for( imp in imp_method) {
          if (imp == 'knn'){
            next
          }
          print(imp)
          imp_method_name <- toString(imp)
          new_data <- integer((numberofzeros * 0.01 * nrow(dataset)))
          temp_Data <- mice(mice_dataframe,m=5,maxit=2,meth=imp,seed=500)
          saved_vector <- integer((numberofzeros * 0.01 * nrow(dataset)))

          # determine the average imputation over all imputed values
          for (p in 1:5) {
            completed_Data <- mice::complete(temp_Data, p)
            new_data <- completed_Data[c(zeros), 1]
            saved_vector <- saved_vector + new_data
          }

          new_data <- saved_vector / 5
          correlation <- cor(saved_data, new_data)

          # Determine if the correlation is below the given cutoff and add the metabolite to the low correlation list if it is
          if (correlation < low_cor_cutoff || is.na(correlation)) {
            na_copy_frame <- mice_dataframe
            na_copy_frame[na_copy_frame == 0] <- NA
            low_cor_vector <- c(low_cor_vector, paste(imp, ": ",zerometa, sep=""))
          }
          df[zerometa, imp] <- round(correlation, 3)
          plotting_df[imp_method_name] <- new_data
        }

        # Here is the KNN method implemented from the VIM package
        if ('knn' %in% imp_method) {
          complete_data2 <- VIM::kNN(mice_dataframe, k=(number_of_correlated_metabolites-1))[,1:number_of_correlated_metabolites]
          new_data2 <- complete_data2[c(zeros), 1]
          correlation <- cor(saved_data, new_data2)
          plotting_df$knn <- new_data2
          if (correlation < low_cor_cutoff || is.na(correlation)) {
            na_copy_frame <- mice_dataframe
            na_copy_frame[na_copy_frame == 0] <- NA
            missingness_count <- length(na.omit(mice_dataframe[,1])) / nrow(dataset) * 100
            print(missingness_count)
            low_cor_vector <- c(low_cor_vector, paste("KNN: ", zerometa, sep=""))
          }
          df[zerometa, "knn"] <- round(cor(saved_data, new_data2),3)  ## add zeros
        }
        # Plotting a percentage of the well performing imputations and a percentage of the poorly imputations.
        data <- melt(plotting_df)
        plotcorrelation <- mean(as.numeric(as.vector(df[zerometa,])), na.rm = TRUE)
        twenty_percent_chance <- sample(1:5, 1)
        if (is.nan(plotcorrelation)){
          plotcorrelation <- 0.0
        }
        if ((plotcorrelation > cutoff) & (twenty_percent_chance==1)) {
          plot(ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.10) + ggtitle(paste("Well performing imputation: ", metabolite, " ", numberofzeros, "%", sep="")))
        }
        fifty_percent_chance <- sample(1:2, 1)
        if ((plotcorrelation < cutoff) & (fifty_percent_chance==1)) {
          plot(ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.10) + ggtitle(paste("Poorly performing imputation: ", metabolite, " ", numberofzeros, "%", sep="")))
        }
      }
    }


    print(paste("Number of variables without enough correlated other variables: ", round(number_of_skipped_metabolites / ncol(dataset) * 100, 2),"%", sep = ""))
    # Check if this list exists and
    if (length(listofskippedmetabolites)) {
      plot(ggplot() + theme_void()+ ggtitle("List of variables without enough correlated other variables."))

      grid.table(listofskippedmetabolites)

    }
    unique_low_cor_vector <- unique(low_cor_vector)

    if (length(unique_low_cor_vector)) {

      plot(ggplot() + theme_void()+ ggtitle("List of variables that performed poorly in the benchmark."))
      grid.table(unique_low_cor_vector)
    }

    # Remove any failed imputations from the dataset.
    df <- na.omit(df)
    # Make a boxplot of the results and identify the best performing method on average.
    colmeans <- colMeans(df[sapply(df, is.numeric)])
    best_method <- names(colmeans[colmeans == max(colmeans)])
    boxplot(df, main=paste("Plot of the correlation of benchmarked methods, with the best method being: ", best_method, sep=""))
    print(paste("The best method for this dataset is: ", best_method, sep=""))
    # Outputs the data into a csv file.
    write.csv(df,'Correlation_per_method.csv')
    dev.off()
    # if setting 1 is selected this is the end of the analysis and the benchmarked correlation will be returned.
    if (impute_setting == 1 ) {
      return(df)
    }
    # Only if setting 2 is selected will the best method be set as the imputation method.
    impute_method <- best_method
  }
  # This part does the imputation on the dataset,
  # on setting 2 the imputation method will be selected from the benchmark, on setting 3 it needs to be supplied
  if (impute_setting == 2 || impute_setting == 3){
    cor_df <- dataset
    cor_df[cor_df == 0] <- NA
    cor_matrix <- as.data.frame(cor(cor_df, method = cor_method))
    metabolite_list <- colnames(dataset)
    number_of_clusters <- ncol(dataset) / 10
    if (ncol(dataset)%%10 != 0){
      number_of_clusters <- number_of_clusters +1
    }

    # This will generate a list of clusters based on correlation in groups of 10 and then removes them from the list
    # untill all variables have been used
    clusters_list <- c()
    for (cluster_number in 1:number_of_clusters) {
      first_metabolite <- metabolite_list[1]
      cor_matrix <- cor_matrix[metabolite_list, ]
      sorteddata <- cor_matrix[order(-cor_matrix[first_metabolite]),]
      cluster_of_correlation <- rownames(sorteddata)[1:10]
      metabolite_list <- setdiff(metabolite_list, cluster_of_correlation)
      clusters_list <- c(clusters_list, cluster_of_correlation)
      clusters_list <- na.omit(unique(clusters_list))
    }

    # Imputes the data using the selected impute method. Starting with the first 10 to make a dataframe.
    mice_dataframe <- (select(cor_df, clusters_list[1:10]))
    if (impute_method == 'knn'){
      completed_df <- VIM::kNN(mice_dataframe, k=9)[, 1:10]
    }
    print(impute_method)
    if (impute_method != 'knn'){
      temp_Data <- mice(mice_dataframe,m=5,maxit=2,meth=impute_method,seed=500)
      completed_df <- mice::complete(temp_Data, 1)
    }

    rownames(completed_df) <- rownames(mice_dataframe)

    # This builds on the previously made dataframe to impute the rest of the dataset.
    oldnum <- 11
    for (newnum in seq(20, length(clusters_list), 10)) {
      print(oldnum:newnum)
      print(newnum)
      print(length(clusters_list))
      if (newnum > length(clusters_list) ){
        next
      }
      clusters <- clusters_list[oldnum:newnum]
      mice_dataframe <- (select(cor_df, clusters))
      if (impute_method == 'knn'){
        completed_Data <- VIM::kNN(mice_dataframe, k=9)[, 1:10]
      }
      print(impute_method)
      if (impute_method != 'knn'){
        temp_Data <- mice(mice_dataframe,m=5,maxit=2,meth=impute_method,seed=500)
        completed_Data <- mice::complete(temp_Data, 1)
      }
      rownames(completed_Data) <- rownames(mice_dataframe)
      completed_df <- merge(completed_df, completed_Data, by=0, all=T)
      rownames(completed_df) <- completed_df$Row.names
      completed_df <- subset(completed_df, select=  -Row.names)
      oldnum <- newnum + 1
    }
    # A correction in case the length of the dataset is not divisible by 10, to still use 10 variables for the imputation
    if (length(clusters_list)%%10 != 0){
      newnum <- length(clusters_list)
      oldnum <- length(clusters_list) - 9
      clusters <- clusters_list[oldnum:newnum]
      mice_dataframe <- (select(cor_df, clusters))
      if (impute_method == 'knn'){
        completed_Data <- VIM::kNN(mice_dataframe, k=9)[, 1:10]
      }
      if (impute_method != 'knn'){
        temp_Data <- mice(mice_dataframe,m=5,maxit=2,meth=impute_method,seed=500)
        completed_Data <- mice::complete(temp_Data, 1)
      }
      rownames(completed_Data) <- rownames(mice_dataframe)
      completed_Data <- completed_Data[, (11-length(clusters_list)%%10):10]
      completed_df <- merge(completed_df, completed_Data, by=0, all=T)
      rownames(completed_df) <- completed_df$Row.names
      completed_df <- subset(completed_df, select=  -Row.names)
    }
    completed_df[completed_df < 0] <- 0.0
    # A couple transformation methods added to give more options for data handling.
    if (transform_method=="none"){
      return(completed_df)
    }
    if (transform_method=="scale"){
      return(apply(completed_df, 2, scale))
    }
    if (transform_method=="minmax"){
      return(apply(completed_df, 2, function (x) { (x - mean(x)) / (max(x) - min(x)) }))
    }
    if (transform_method=="log"){
      return(apply(completed_df, 2,log))
    }
    if (transform_method=="zscore"){
      return(apply(completed_df, 2, function (x) { (x - mean(x)) / stats::sd(x) }))
    }
  }
}
