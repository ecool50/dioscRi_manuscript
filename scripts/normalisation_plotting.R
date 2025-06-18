# Return a matrix of values (mean or maybe var).
# Rows are cytokines (or any marker).
# Cols are anchor samples.
# Does it matter if trans=FALSE or TRUE?
# Similar to above, but not subpopulations. All cell events.
get_summaries_per_channel <- function(data, markers, samples){
    outmat <- c()
    for(sample in samples){
        data_sub <- data %>%
            dplyr::filter(sample_id == sample)
        vals <- apply(data_sub[,markers], MARGIN=2, FUN="mean"); # or var
        outmat <- cbind(outmat, vals);
        colnames(outmat)[ncol(outmat)] <- sample;
    }
    return(outmat);
}

get_summaries_per_celltype <- function(data, markers, celltypes){
    outmat <- c()
    for(celltype in celltypes){
        data_sub <- data %>%
            dplyr::filter(Celltype == celltype)
        vals <- apply(data_sub[,markers], MARGIN=2, FUN="mean"); # or var
        outmat <- cbind(outmat, vals);
        colnames(outmat)[ncol(outmat)] <- celltype;
    }
    return(outmat);
}

# Compute the trace of the covariance matrix.
# Input matricies are  Rows:markers   X   Cols:barcode samples
varDiffTotalPrePostT <- function(tbl_pre, tbl_post){
    
    ####   Transpose   ####
    tbl_pre <- t(tbl_pre);
    tbl_post <- t(tbl_post);
    #######################
    
    cpre <- cov(tbl_pre);
    cpost <- cov(tbl_post);
    epre <- eigen(cpre, symmetric=TRUE, only.values=TRUE);
    epost <- eigen(cpost, symmetric=TRUE, only.values=TRUE);
    trace_pre <- sum(epre$values);
    trace_post <- sum(epost$values);
    return(trace_pre - trace_post);
}

# Swap columns of pre/post: 1to1, 2to2, etc
# Calculate varDiffTotal for all possible permutations. Return vector.
permuteColsPrePost <- function(mat_pre, mat_post){
    if(!(all(colnames(mat_pre) == colnames(mat_post)))){
        stop("pre/post colnames don't match");
    }
    mat <- cbind(mat_pre, mat_post);
    N_bcs <- ncol(mat_pre);
    print(N_bcs)
    tfvec <- rep(1, N_bcs);
    count <- 0;
    testres <- rep(1, (2^N_bcs));
    for(itr in 0:((2^N_bcs)-1)){
        tfvec1 <- as.logical(intToBits(itr)[1:N_bcs]);
        tfvec <- c(tfvec1, !tfvec1)
        perm_mat_pre <- mat[,tfvec];
        perm_mat_post <- mat[,!tfvec];
        count <- count + 1;
        testres[count] <- varDiffTotalPrePostT(perm_mat_pre, perm_mat_post);
    }
    return(testres);
}

library(parallel)

permuteColsPrePost_parallel <- function(mat_pre, mat_post) {
    if (!(all(colnames(mat_pre) == colnames(mat_post)))) {
        stop("pre/post colnames don't match")
    }
    
    mat <- cbind(mat_pre, mat_post)
    N_bcs <- ncol(mat_pre)
    print(N_bcs)
    
    # Total number of iterations
    total_permutations <- 2^N_bcs
    
    # Define the task for each worker
    compute_permutation <- function(itr) {
        tfvec1 <- as.logical(intToBits(itr)[1:N_bcs])
        tfvec <- c(tfvec1, !tfvec1)
        perm_mat_pre <- mat[, tfvec]
        perm_mat_post <- mat[, !tfvec]
        varDiffTotalPrePostT(perm_mat_pre, perm_mat_post)
    }
    
    # Set up parallel processing
    cl <- makeCluster(detectCores() - 1) # Use one less than the total cores
    clusterExport(cl, varlist = c("mat", "N_bcs", "varDiffTotalPrePostT"), envir = environment())
    clusterEvalQ(cl, library(parallel)) # Ensure the library is loaded on each worker
    
    # Run the computation in parallel
    testres <- parSapply(cl, 0:(total_permutations - 1), compute_permutation)
    
    stopCluster(cl) # Stop the cluster
    
    return(testres)
}


# boxplot (top left)
# rows: cytokines, cols: anchor samples
# addPoints: add a point to the boxplot for each sample. (may want to adjust point size)
boxSummaryValsPrePost <- function(mat_pre, mat_post, colorPre="blue", 
                                  colorPost="wheat", addPoints=FALSE){
    # make a list
    valueslist <- list();
    for(rowi in nrow(mat_pre):1){
        vname <- sprintf("post %s", rownames(mat_post)[rowi]);
        valueslist[[vname]] <- mat_post[rowi,];
        vname <- sprintf("pre %s", rownames(mat_pre)[rowi]);
        valueslist[[vname]] <- mat_pre[rowi,];
    }
    boxplot(valueslist, boxwex=.5, xaxt="n", yaxt="n", cex.axis=1, las=2, 
            range=0, lwd=1, boxlwd=1.25, horizontal=TRUE, 
            col=rep(c(colorPost, colorPre), 2*nrow(mat_pre)), labels=T);
    axis(side=1, las=0); # 1=below
    
    if(addPoints){
        if(nrow(mat_pre) > 15){
            ptcex <- .3;
        } else{
            ptcex <- .8;
        }
        for(vi in 1:length(valueslist)){
            points(y=rep(vi, length(valueslist[[vi]])), x=valueslist[[vi]], 
                   cex=ptcex, lwd=1, col="darkred");
            points(y=rep(vi, length(valueslist[[vi]])), x=valueslist[[vi]], 
                   cex=ptcex/2, lwd=1, col="yellow");
        }
    }
}


# boxplot (top left)
# rows: cytokines, cols: anchor samples
# addPoints: add a point to the boxplot for each sample. (may want to adjust point size)
boxSummaryValsPrePostCelltype <- function(mat_pre, mat_post, colorPre="blue", 
                                  colorPost="wheat", addPoints=FALSE){
    
    mat_pre <- t(mat_pre)
    mat_post <- t(mat_post)
    # make a list
    valueslist <- list();
    for(rowi in nrow(mat_pre):1){
        vname <- sprintf("post %s", rownames(mat_post)[rowi]);
        valueslist[[vname]] <- mat_post[rowi,];
        vname <- sprintf("pre %s", rownames(mat_pre)[rowi]);
        valueslist[[vname]] <- mat_pre[rowi,];
    }
    boxplot(valueslist, boxwex=.5, xaxt="n", yaxt="n", cex.axis=1, las=2, 
            range=0, lwd=1, boxlwd=1.25, horizontal=TRUE, 
            col=rep(c(colorPost, colorPre), 2*nrow(mat_pre)), labels=T);
    axis(side=1, las=0); # 1=below
    
    if(addPoints){
        if(nrow(mat_pre) > 15){
            ptcex <- .3;
        } else{
            ptcex <- .8;
        }
        for(vi in 1:length(valueslist)){
            points(y=rep(vi, length(valueslist[[vi]])), x=valueslist[[vi]], 
                   cex=ptcex, lwd=1, col="darkred");
            points(y=rep(vi, length(valueslist[[vi]])), x=valueslist[[vi]], 
                   cex=ptcex/2, lwd=1, col="yellow");
        }
    }
}


varBarPlot <- function(mat_pre, mat_post, colorPre="blue", colorPost="wheat"){
    if(nrow(mat_pre) > 15){
        cexnames <- .7;
    } else{
        cexnames <- 1;
    }
    vpre <- apply(mat_pre, MARGIN=1, FUN=var);
    vpost <- apply(mat_post, MARGIN=1, FUN=var);
    mat <- rbind( rev(vpost), rev(vpre));
    barplot(height=mat, horiz=TRUE, beside=TRUE, col=rep(c(colorPost, colorPre), 2*nrow(mat_pre)), las=1, cex.names=cexnames);
    #title(main="Variance (in mean signal)");
    title(main="Variance");
}

varBarPlotCelltype <- function(mat_pre, mat_post, colorPre="blue", colorPost="wheat"){
    if(nrow(mat_pre) > 15){
        cexnames <- .7;
    } else{
        cexnames <- 1;
    }
    vpre <- apply(t(mat_pre), MARGIN=1, FUN=var);
    vpost <- apply(t(mat_post), MARGIN=1, FUN=var);
    mat <- rbind( rev(vpost), rev(vpre));
    barplot(height=mat, horiz=TRUE, beside=TRUE, col=rep(c(colorPost, colorPre), 2*nrow(mat_pre)), las=1, cex.names=cexnames);
    #title(main="Variance (in mean signal)");
    title(main="Variance");
}

# Three plots: 
#  1. horizontal boxplot, pre and post box for each channel.
#  2. Variance of plot 1. 
#  3. Null distribution and p-value
# Total variance (reduction) in mean cytokine levels for all cell events
# Permutation test for significance swapping pre/post.
totalVar_allEvents <- function(mat_pre, mat_post, markers, colorPre="blue", colorPost="wheat"){
    cols_to_use <- markers
    if(!(all(colnames(mat_pre) == colnames(mat_post)))){
        stop("totalVar_allEvents: File names in pre and post directories don't match.");
    }
    test_real <- varDiffTotalPrePostT(mat_pre, mat_post);
    perm_tests <- permuteColsPrePost(mat_pre, mat_post);
    N_as_or_more_extreme <- sum(perm_tests >= test_real);
    pv <- N_as_or_more_extreme / length(perm_tests);
    
    pngwidth<-1600; pngheight<-1600;
    png(filename="../Plots/PrePostVarianceAllCelltype.png", width=pngwidth, height=pngheight);
    layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE));
    
    title <- sprintf("%s", "Mean signal intensity");
    
    # plot 1:
    par(cex=1.5);
    par(mar=c(2,1,2,0)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
    boxSummaryValsPrePost(mat_pre, mat_post, colorPre, colorPost);
    title(main=title);
    legend("bottomright", legend=c("Pre", "Post"), col=c(colorPre, colorPost), lwd=9);
    legend("bottomright", legend=c("Pre", "Post"), col=c(colorPre, colorPost), lwd=9);
    
    # plot 2: variance bars
    #par(mar=c(2,0,2,1)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
    par(mar=c(2,7.5,2,1)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
    varBarPlot(mat_pre, mat_post, colorPre, colorPost);
    
    # plot 3: Bottom NULL dist hist
    par(cex=1.5);
    par(mar=c(5,4.5,4,1)+0.1); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
    #hist(perm_tests, breaks=70, main="", xlab="Test statistic null distribution");
    #hist(perm_tests, breaks=70, main="", xlab="Change in total variance null distribution");
    hist(perm_tests, breaks=70, main="", xlab="Change in total variance", cex.lab=1.5);
    abline(v=test_real, col=2, lwd=5);
    #title(main=sprintf("p = %.05f", pv), line=-1);
    title(main=sprintf("p = %.05f", pv), line=0);
    dev.off();
    print(sprintf("p = %.05f", pv), q=F);
}

# Three plots: 
#  1. horizontal boxplot, pre and post box for each channel.
#  2. Variance of plot 1. 
#  3. Null distribution and p-value
# Total variance (reduction) in mean cytokine levels for all cell events
# Permutation test for significance swapping pre/post.
totalVar_allEventsCelltype <- function(mat_pre, mat_post, celltypes, colorPre="blue", colorPost="wheat"){
    
    cols_to_use <- celltypes
    # mat_pre <- t(mat_pre)
    # mat_post <- t(mat_post)
    if(!(all(colnames(mat_pre) == colnames(mat_post)))){
        stop("totalVar_allEvents: File names in pre and post directories don't match.");
    }
    test_real <- varDiffTotalPrePostT(t(mat_pre), t(mat_post));
    perm_tests <- permuteColsPrePost_parallel(t(mat_pre), t(mat_post));
    N_as_or_more_extreme <- sum(perm_tests >= test_real);
    pv <- N_as_or_more_extreme / length(perm_tests);
    
    pngwidth<-1600; pngheight<-1600;
    png(filename="../Plots/PrePostVarianceAllMGcelltypes.png", width=pngwidth, height=pngheight);
    layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE));
    
    title <- sprintf("%s", "Combined - mean signal intensity per cell type");
    
    # plot 1:
    par(cex=1.5);
    par(mar=c(2,1,2,0)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
    boxSummaryValsPrePostCelltype(mat_pre, mat_post, colorPre, colorPost);
    title(main=title);
    legend("bottomright", legend=c("Pre", "Post"), col=c(colorPre, colorPost), lwd=9);
    legend("bottomright", legend=c("Pre", "Post"), col=c(colorPre, colorPost), lwd=9);
    
    # plot 2: variance bars
    #par(mar=c(2,0,2,1)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
    par(mar=c(2,7.5,2,1)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
    varBarPlotCelltype(mat_pre, mat_post, colorPre, colorPost);
    
    # plot 3: Bottom NULL dist hist
    par(cex=1.5);
    par(mar=c(5,4.5,4,1)+0.1); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
    #hist(perm_tests, breaks=70, main="", xlab="Test statistic null distribution");
    #hist(perm_tests, breaks=70, main="", xlab="Change in total variance null distribution");
    hist(perm_tests, breaks=70, main="", xlab="Change in total variance", cex.lab=1.5);
    abline(v=test_real, col=2, lwd=5);
    #title(main=sprintf("p = %.05f", pv), line=-1);
    title(main=sprintf("p = %.05f", pv), line=0);
    dev.off();
    print(sprintf("p = %.05f", pv), q=F);
}

boxSummaryValsPrePost_ggplot <- function(raw, dioscri, cycombine, imubac, type = "both",
                                         legendPos = "none", title = "", addPoints=FALSE) {
    # Create a combined dataframe for plotting
    raw_df <- as.data.frame(t(raw))
    dioscRi_df <- as.data.frame(t(dioscri))
    cycombine_df <- as.data.frame(t(cycombine))
    imubac_df <- as.data.frame(t(imubac))
    # scmerge2_df <- as.data.frame(t(scmerge2))
    
    # Melt the dataframes to long format for ggplot
    raw_long <- reshape2::melt(raw_df, variable.name = "Marker", value.name = "Value")
    dioscRi_long <- reshape2::melt(dioscRi_df, variable.name = "Marker", value.name = "Value")
    cycombine_long <- reshape2::melt(cycombine_df, variable.name = "Marker", value.name = "Value")
    imubac_long <- reshape2::melt(imubac_df, variable.name = "Marker", value.name = "Value")
    # scmerge2_long <- melt(scmerge2_df, variable.name = "Marker", value.name = "Value")
    
    # Add labels for raw and dioscRi
    raw_long$Data <- "Raw"
    dioscRi_long$Data <- "dioscRi"
    cycombine_long$Data <- "cyCombine"
    imubac_long$Data <- "iMUBAC"
    # scmerge2_long$Data <- "scMerge2"
    
    # Combine raw and dioscRi data
    combined_df <- rbind(raw_long, dioscRi_long, cycombine_long, 
                         imubac_long) %>%
        dplyr::mutate(Data = factor(Data, 
                                         levels = c("dioscRi","cyCombine","iMUBAC","Raw")))
    
    if(type == "both"){
        # Create the plot
        plot <- combined_df %>%
            dplyr::filter(Data %in% c("dioscRi", "Raw")) %>% 
            ggplot(aes(x = Marker, y = Value, fill = Data)) +
            geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +
            # scale_fill_manual(values = c("Raw" = colorraw, "dioscRi" = colordioscRi)) +
            coord_flip() + # Flip the axes for horizontal boxplots
            theme_bw() +
            ggtitle(title) +
            theme(
                axis.text.y = element_text(size = 12),
                axis.text.x = element_text(size = 12),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12),
                legend.position = legendPos,
                legend.title = element_blank(),
                xis.line = element_line(colour = "black"),
                # panel.grid.major = element_blank(),
                # panel.grid.minor = element_blank(),
                # panel.border = element_blank(),
                # panel.background = element_blank(),
                plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
            ) +
            labs(x = "", y = "") +
            scale_fill_npg()
    } else{
        # Create the plot
        plot <- combined_df %>%
            ggplot(aes(x = Marker, y = Value, fill = Data)) +
            geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +
            # scale_fill_manual(values = c("Raw" = colorraw, "dioscRi" = colordioscRi)) +
            coord_flip() + # Flip the axes for horizontal boxplots
            theme_bw() +
            ggtitle(title) +
            theme(
                axis.text.y = element_text(size = 12),
                axis.text.x = element_text(size = 12),
                axis.title.x = element_text(size = 12),
                axis.title.y = element_text(size = 12),
                legend.position = legendPos,
                legend.title = element_blank(),
                # xis.line = element_line(colour = "black"),
                # panel.grid.major = element_blank(),
                # panel.grid.minor = element_blank(),
                # panel.border = element_blank(),
                # panel.background = element_blank(),
                plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
            ) +
            labs(x = "", y = "") +
            scale_fill_npg()
    }
    
    
    
    
    return(plot)
}

varBarPlot_ggplot <- function(raw, dioscri, cycombine, imubac,
                              legendPos = "none", title = "", type = "both",
                              colorPre = "wheat", colorPost = "blue") {
    # Calculate variance for each row in pre and post matrices
    vraw <- apply(raw, 1, var)
    vdioscri <- apply(dioscri, 1, var)
    vcycombine <- apply(cycombine, 1, var)
    vimubac <- apply(imubac, 1, var)
    # vscmerge2 <- apply(scmerge2, 1, var)
    
    # Create a data frame for ggplot
    variance_df <- data.frame(
        Row = factor(rownames(raw), levels = rownames(raw)), # Reverse row order
        Variance = c(vraw, vdioscri, vcycombine, vimubac),
        Condition = rep(c("Raw", "dioscRi", "cyCombine", "iMUBAC"), 
                        each = nrow(raw))
    ) %>%
        dplyr::mutate(Condition = factor(Condition, 
                                         levels = c("dioscRi","cyCombine","iMUBAC","Raw")))
    
    # Create the bar plot
    if(type == "both"){
        plot <- variance_df %>%
            dplyr::filter(Condition %in% c("dioscRi", "Raw")) %>%
        ggplot(aes(x = Variance, y = Row, fill = Condition),
               position = position_dodge(width = 0.9)) +
            geom_bar(stat = "identity", position = "dodge") +
            # scale_fill_manual(values = c("Raw" = colorPre, "dioscRi" = colorPost)) +
            theme_bw() +
            theme(
                axis.text.y = element_text(size = 12), # Adjust text size
                axis.text.x = element_text(size = 12),
                axis.title = element_text(size = 12),
                legend.position = legendPos,
                legend.title = element_blank(),
                xis.line = element_line(colour = "black"),
                # panel.grid.major = element_blank(),
                # panel.grid.minor = element_blank(),
                # panel.border = element_blank(),
                # panel.background = element_blank(),
                plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
            ) +
            labs(title = title, x = "", y = "") +
            scale_fill_npg()
    } else{
        plot <- variance_df %>%
            ggplot(aes(x = Variance, y = Row, fill = Condition),
                   position = position_dodge(width = 0.9)) +
            geom_bar(stat = "identity", position = "dodge") +
            # scale_fill_manual(values = c("Raw" = colorPre, "dioscRi" = colorPost)) +
            theme_bw() +
            theme(
                axis.text.y = element_text(size = 12), # Adjust text size
                axis.text.x = element_text(size = 12),
                axis.title = element_text(size = 12),
                legend.position = legendPos,
                legend.title = element_blank(),
                xis.line = element_line(colour = "black"),
                # panel.grid.major = element_blank(),
                # panel.grid.minor = element_blank(),
                # panel.border = element_blank(),
                # panel.background = element_blank(),
                plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
            ) +
            labs(title = title, x = "", y = "") +
            scale_fill_npg()
    }
    
    
    return(plot)
}

cvBoxPlotByMarker_ggplot <- function(raw, dioscri, cycombine, imubac,
                                     legend_pos = "right",
                                     title = "Relative CV per Marker") {
    # Helper to melt and label each matrix
    long_data <- function(mat, label) {
        if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
        
        df <- as.data.frame(mat)
        df$marker <- rownames(df)
        
        reshape2::melt(df, id.vars = "marker",
                       variable.name = "sample",
                       value.name = "value") |>
            dplyr::mutate(method = label)
    }
    
    # Combine long-form data from all methods
    df_all <- dplyr::bind_rows(
        long_data(raw, "Raw"),
        long_data(dioscri, "dioscRi"),
        long_data(cycombine, "cyCombine"),
        long_data(imubac, "iMUBAC")
    )
    
    # Compute CV per marker per method
    cv_df <- df_all |>
        dplyr::group_by(marker, method) |>
        dplyr::summarise(
            mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            cv = ifelse(mean != 0, sd / mean, NA_real_),
            .groups = "drop"
        )
    
    # Extract raw CVs for reference
    raw_cv <- cv_df |>
        dplyr::filter(method == "Raw") |>
        dplyr::select(marker, raw_cv = cv)
    
    # Join and compute relative CVs
    cv_df <- cv_df |>
        dplyr::left_join(raw_cv, by = "marker") |>
        dplyr::mutate(relative_cv = cv / raw_cv) |>
        dplyr::filter(method != "Raw")  # Optional: drop Raw if you only want relatives
    
    # Set factor levels for consistent method order
    cv_df$method <- factor(cv_df$method,
                           levels = c("dioscRi", "cyCombine", "iMUBAC"))
    
    # Plot relative CVs
    ggplot(cv_df, aes(x = method, y = relative_cv, fill = method)) +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
        # coord_flip() +
        labs(
            title = title,
            x = "",
            y = "Relative CV to Raw"
        ) +
        theme_bw() +
        theme(
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.position = "none",
            legend.title = element_blank(),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            aspect.ratio = 1
        ) + 
        scale_fill_npg() 
}

#' Plot multiple ROC curves and a separate AUC barplot
#'
#' @param rocs A named list of 'roc' objects from the pROC package.
#' @param breaks Tick marks on axes.
#' @param legend_title Title for the legend.
#' @param data_name Title for the ROC plot.
#' @param legend_order Optional character vector for legend and bar order.
#'
#' @return A list with `roc_plot` and `auc_bar`.
ggrocs <- function(
        rocs,
        breaks       = seq(0, 1, 0.1),
        legend_title = "Legend",
        data_name    = "",
        legend_order = NULL,
        color_scale  = NULL,
        which_plot   = c("both", "roc", "bar")
) {
    # —— argument checking
    which_plot <- match.arg(which_plot)
    if (length(rocs) == 0) stop("No ROC objects provided.")
    if (!requireNamespace("plyr", quietly = TRUE)) {
        stop("The 'plyr' package is required.")
    }
    library(ggplot2); library(plyr)
    
    # —— compute AUCs & names
    aucs         <- sapply(rocs, function(x) x$auc)
    method_names <- names(rocs)
    if (!is.null(legend_order)) {
        method_names <- intersect(legend_order, method_names)
    }
    
    # —— pretty labels only for ROC legend
    pretty_labels <- sprintf("%s (AUC = %.2f)", method_names, aucs[method_names])
    label_map     <- setNames(pretty_labels, method_names)
    
    # —— build ROC data.frame
    roc_vals <- plyr::ldply(method_names, function(m) {
        roc_obj <- rocs[[m]]
        if (!inherits(roc_obj, "roc")) {
            stop("All elements must be 'roc' objects from the pROC package.")
        }
        data.frame(
            fpr    = rev(roc_obj$specificities),
            tpr    = rev(roc_obj$sensitivities),
            method = m,
            stringsAsFactors = FALSE
        )
    })
    roc_vals$method <- factor(roc_vals$method, levels = method_names)
    roc_vals$label  <- factor(label_map[as.character(roc_vals$method)],
                              levels = pretty_labels)
    
    # —— build AUC bar data.frame (raw names, no AUC in factor)
    auc_df <- data.frame(
        Method = factor(method_names, levels = method_names),
        AUC    = as.numeric(aucs[method_names])
    )
    
    # —— optional: split colour_scale into raw vs. pretty
    if (!is.null(color_scale)) {
        raw_colors    <- color_scale[method_names]
        pretty_colors <- setNames(raw_colors, pretty_labels)
    }
    
    # —— ROC plot
    roc_plot <- ggplot(roc_vals, aes(x = fpr, y = tpr, colour = label)) +
        geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0),
                     colour = "grey80", linetype = "dashed", size = 1) +
        geom_line(size = 1) +
        scale_x_reverse(
            name   = "False Positive Rate (1 - Specificity)",
            limits = c(1, 0), breaks = breaks
        ) +
        scale_y_continuous(
            name   = "True Positive Rate (Sensitivity)",
            limits = c(0, 1), breaks = breaks
        ) +
        theme(aspect.ratio = 1, plot.margin = margin(2, 2, 2, 2)) +
        ggtitle(data_name) +
        guides(colour = guide_legend(legend_title)) +
        theme(
            plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5),
            axis.title.x = element_text(color = "black", size = 12),
            axis.title.y = element_text(color = "black", size = 12),
            axis.text = element_text(size = 12),
            # legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.position = c(1, 0),
            legend.justification = c(1, 0),
            legend.background = element_rect(fill = "transparent"),
            panel.background = element_blank(),
            axis.ticks = element_line(color = "grey80")
        )
    
    if (!is.null(color_scale)) {
        roc_plot <- roc_plot +
            scale_colour_manual(
                values = pretty_colors,
                breaks = pretty_labels    # enforce order
            )
    }
    
    # —— AUC bar plot
    auc_bar <- ggplot(auc_df, aes(x = Method, y = AUC, fill = Method)) +
        geom_col(width = 0.3) +
        geom_text(aes(label = round(AUC, 2)), vjust = -0.5, size = 4) +
        scale_y_continuous(expand = c(0, 0.1)) +
        labs(y = "AUC", x = NULL, title = data_name) +
        theme_minimal(base_size = 12) +
        theme(
            plot.title     = element_text(face = "bold", hjust = 0.5),
            axis.text.x    = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            aspect.ratio    = 1
        )
    
    if (!is.null(color_scale)) {
        auc_bar <- auc_bar +
            scale_fill_manual(
                values = raw_colors,
                breaks = method_names
            )
    }
    
    # —— return
    switch(
        which_plot,
        roc = roc_plot,
        bar = auc_bar,
        both = list(roc_plot = roc_plot, auc_bar = auc_bar)
    )
}

mse_barplot_by_method <- function(raw, dioscri, cycombine, imubac,
                                  title = "Mean Squared Error (vs. Raw)") {
    # Helper: melt matrix and add label
    long_data <- function(mat, label) {
        if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
        
        df <- as.data.frame(mat)
        df$marker <- rownames(df)
        reshape2::melt(df, id.vars = "marker", variable.name = "sample", value.name = "value") |>
            dplyr::mutate(method = label)
    }
    
    # Melt all matrices
    df_raw <- long_data(raw, "Raw")
    df_dioscri <- long_data(dioscri, "dioscRi")
    df_cycombine <- long_data(cycombine, "cyCombine")
    df_imubac <- long_data(imubac, "iMUBAC")
    
    # Combine all and align on marker+sample
    all_methods <- dplyr::bind_rows(df_dioscri, df_cycombine, df_imubac)
    reference <- df_raw |>
        dplyr::select(marker, sample, raw_value = value)
    
    mse_df <- all_methods |>
        dplyr::left_join(reference, by = c("marker", "sample")) |>
        dplyr::mutate(sq_error = (value - raw_value)^2) |>
        dplyr::group_by(method) |>
        dplyr::summarise(
            mse = mean(sq_error, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Plot
    ggplot(mse_df, aes(x = method, y = mse, fill = method)) +
        geom_col(width = 0.6) +
        geom_text(aes(label = round(mse, 4)), vjust = -0.5, size = 4) +
        labs(x = "", y = "Mean Squared Error", title = title) +
        scale_fill_npg() +
        theme_bw() +
        theme(
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.position = "none",
            aspect.ratio = 1
        )
}

mse_boxplot_by_marker <- function(raw, dioscri, cycombine, imubac,
                                  title = "Mean Squared Error per Marker (vs. Raw)") {
    library(reshape2)
    library(dplyr)
    library(ggplot2)
    library(ggsci)
    
    # Melt helper with synthetic cell_id
    melt_matrix <- function(mat, method_label) {
        if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
        df <- as.data.frame(mat)
        df$cell_id <- seq_len(nrow(df))  # use row number as cell_id
        melt(df, id.vars = "cell_id", variable.name = "marker", value.name = "value") |>
            mutate(method = method_label)
    }
    
    # Melt all matrices
    df_raw       <- melt_matrix(raw, "Raw")
    df_dioscri   <- melt_matrix(dioscri, "dioscRi")
    df_cycombine <- melt_matrix(cycombine, "cyCombine")
    df_imubac    <- melt_matrix(imubac, "iMUBAC")
    
    # Stop if raw is invalid
    if (is.null(df_raw)) stop("Raw matrix is required.")
    
    # Combine all methods (non-null)
    all_methods <- bind_rows(df_dioscri, df_cycombine, df_imubac)
    if (nrow(all_methods) == 0) stop("No valid method matrices provided.")
    
    # Raw values for matching
    df_ref <- df_raw |> 
        dplyr::select(cell_id, marker, raw_value = value)
    
    # Join and compute MSE per marker
    mse_df <- all_methods |>
        dplyr::left_join(df_ref, by = c("cell_id", "marker")) |>
        dplyr::mutate(sq_error = (value - raw_value)^2) |>
        dplyr::group_by(marker, method) |>
        dplyr::summarise(mse = mean(sq_error, na.rm = TRUE), .groups = "drop")
    
    # Optional: enforce method order
    method_levels <- c("dioscRi", "cyCombine", "iMUBAC")
    mse_df$method <- factor(mse_df$method, levels = method_levels[method_levels %in% mse_df$method])
    
    # Compute mean MSE per method
    mean_labels <- mse_df %>%
        group_by(method) %>%
        summarise(mean_mse = mean(mse, na.rm = TRUE), .groups = "drop")
    
    print(mean_labels)
    
    
    # Plot MSE as boxplot across markers
    ggplot(mse_df, aes(x = method, y = mse, fill = method)) +
        geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(
            title = title,
            x = "",
            y = "Mean Squared Error vs. Raw"
        ) +
        scale_fill_npg() +
        theme_bw() +
        theme(
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.position = "none",
            aspect.ratio = 1
        )
}