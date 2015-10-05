#' plotOutliers
#'
#' This function plot the outliers and non-outliers between using the table variables
#'
#' @param  file the file that contains the expeted value and predicted values
#'

plotOutliers <- function(protFile, height = 800, width = 800, cols = 3){

    # Read the data files.
    proData <- read.table(file=protFile,  header = TRUE, sep = "\t")

    # Remove the empty values and other data, also remove the string data Proteins and Peptides
    datPR <- removeFirstColumn(proData)

    datPR <- processData(datPR)

    plotsOut <- plotOutlierDistribution(datPR)

    png("outliersDistributionPlots.png", width = 800, height = 800)
    multiplot(plotlist = plotsOut, cols=2)
    dev.off()

    plotsOut1 <- plotOutlierOverall(datPR)

    png("outliersOverallPlots.png", width = 800, height = 800)
    multiplot(plotlist = plotsOut1, cols=2)
    dev.off()

}


#' plotOutlierDistribution
#'
#' This function plot the outliers and non-outliers population distribution between using the table variables
#'
#' @param  data the data.frame that contains the expeted value and predicted values
#'

plotOutlierDistribution <- function(data, na.rm = TRUE, ...) {

    #getting subset to analysis, only a few of predictors
    data <- subset.data.frame(data, select = c(EXPERIMENTAL, ITERATIVE_LEHNINGER, BJELL_EXPASY, COFACTOR, SVM))

    nm <- names(data)
    #    breaks <- seq(from = 0, to = 14, by = 0.5)
    #    print(breaks)
    plots <- list()  # new empty list
    for (i in nm)
    {
        if(data[1L] != data[i]){

            newData <- data.frame(x = data[1L], y = data[i])
            colnames(newData) <- c("x", "y")
            pData <- mutate(newData, outlier = (abs(x - y) <= 2*sd(y))) #criteria to detect outliers(if condition is TRUE: non-outlier, else if FALSE: outlier)

            pData[pData=="TRUE"]  <- "non-outliers"    #rename variable to more descriptive name
            pData[pData=="FALSE"] <- "outliers"


            #computing percent outliers: %outlier = FALSE*100/(TRUE + FALSE)
            f_outliers      <-  factor(pData$outlier)
            freq_outliers   <-  count(f_outliers) #retrieve frecuency table into dataFrame
            freq_outliers   <-  mutate(freq_outliers, percent= freq/sum(freq)*100) #add % outliers values
            percent_value   <-  freq_outliers$percent[2] #getting percent of outliers(FALSE)

            plot <- ggplot(pData, aes(x=y, y=..density.., fill=outlier)) +

                geom_histogram(position="identity", alpha=0.3) +
                #geom_density(position="identity", alpha=0.3) +
                geom_line(stat="density", position = "identity", linetype = 2)+
                scale_fill_manual(values = c("black","light gray"))+
                #facet_grid(outlier ~ .,scales="free")+
                labs(fill="Distribution")+
                ggtitle(colnames(data[i]))+
                annotate("text", label=parseOutliersValue(percent_value), parse=TRUE, x=-Inf, y=Inf, hjust=-.2, vjust=2)+
                theme(plot.title=element_text(size=rel(0.6), lineheight=.6))+
                theme_bw()+
                xlab("Isoelectric Point")+
                theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
            plots[[i]] <- plot  # add each plot into plot list
        }

    }
    return (plots)
}


#' plotOutlierOverall
#'
#' This function plot the outliers and non-outliers dispersion between using the table variables
#'
#' @param  data the data.frame that contains the expeted value and predicted values
#'

plotOutlierOverall <- function(data, na.rm = TRUE, ...) {

    #getting subset to analysis
    data <- subset.data.frame(data, select = c(EXPERIMENTAL, ITERATIVE_LEHNINGER, BJELL_EXPASY, COFACTOR, SVM))
    nm <- names(data)

    plots <- list()  # new empty list
    for (i in nm)
    {
        if(data[1L] != data[i]){

            newData <- data.frame(x = data[1L], y = data[i])
            colnames(newData) <- c("x", "y")
            pData <- mutate(newData, outlier = (abs(x - y) <= 2*sd(y))) #criteria to detect outliers

            pData[pData=="TRUE"]  <- "non-outliers"    #rename variable to more descriptive name
            pData[pData=="FALSE"] <- "outliers"


            #computing % outliers
            f_outliers      <-  factor(pData$outlier)
            freq_outliers   <-  count(f_outliers) #retrieve frecuency table into dataFrame
            freq_outliers   <-  mutate(freq_outliers, percent= freq/sum(freq)*100) #add % outliers values
            percent_value   <-  freq_outliers$percent[2] #getting percent of outliers


            plot <- ggplot(pData, aes(x=x, y=y, shape=outlier))+
                geom_point(size=2.0, alpha=.4)+                           # Set points size and transparency
                scale_shape_manual(values=c(1,4)) +                       # Set points shape
                scale_colour_hue(l=50) +                                  # Use a slightly darker palette than normal
                geom_smooth(aes(group=1), method=lm, se=FALSE) +          # Add shaded confidence region
                xlab("EXPERIMENTAL") +
                ylab(i) +
                #facet_grid(outlier ~ .,scales="free")+
                annotate("text", label=parseOutliersValue(percent_value), parse=TRUE, x=-Inf, y=Inf, hjust=-.2, vjust=2)+
                theme_bw() +
                guides(shape=guide_legend(title=NULL))+
                theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
            plots[[i]] <- plot  # add each plot into plot list
        }

    }
    return (plots)
}




#' parseOutliersValue
#'
#' This function parse numerical value into formule or expression
#'
#' @param  value the value to include into any expression
#'

parseOutliersValue <- function(value){
    sy <- as.character("% outliers")
    value = format(value, digits=3)
    eqn <- as.character(as.expression(
        substitute(italic('outliers %') == out,
                   list(out = value))))
    return (eqn)
}
