#' plotPTMShiftOnFractions
#'
#' This fucntion plots the shift on isoelectric point estimation produt of any modification (PTM) of all the predicted varaibles (theoretical mean of the fraction) vs the expted variable (first column). It plots error bars.
#' Both data1 as data2 must have the same size and shape.
#'
#' @param data1 the data.frame with estimations excluding PTM contribution
#' @param data2 the data.frame with estimations including PTM contribution
#'
#'

plotPTMShiftOnFractions <- function (dat1, dat2){
    nm <- names(dat1)
    plots <- list()  # new empty list
    for (i in nm)
    {
        if(dat1[1L] != dat1[i])
        {
            newData <- data.frame(x = dat1[1L], y1 = dat1[i], y2 = dat2[i])
            colnames(newData) <- c("x", "y1", "y2") #x: x-Axis, y1: y-Axis from file1 and y2: y-Axis from file2
            fractions_stats <- ddply(newData, ~x, summarise, y_mean1=mean(y1), y_sd1=sd(y1), y_mean2=mean(y2), y_sd2=sd(y2)) # split data in fraction. compute mean and standard deviation
            dataStats <- reshape(fractions_stats, idvar= "x", varying= list(c("y_mean1","y_mean2"),c("y_sd1","y_sd2")), #reshape data to plot into groups (i.e. modified and non-modified)
                                 v.names = c("y","sd"), times=c("excluding modification","including modification"), direction = "long")
            #print(fractions_stats)
            plot <- ggplot(dataStats, aes(x=x, y=y, shape=time))+
                geom_errorbar(aes(ymin=y-sd, ymax= y+sd), width=.1, size=0.2, colour="black") +
                geom_point(size=2.5, alpha=.4)+                          #set points size and transparency
                scale_shape_manual(values=c(1,2)) +                      #set points shape
                scale_y_continuous(breaks=seq(4, 10, 1)) +
                scale_x_continuous(breaks=seq(4, 7, .25))+
                xlab (nm[1L]) +
                ylab(i) +
                theme_bw()+
                theme(legend.position="none") +
                #theme(legend.background=element_blank()) +
                #theme(legend.key=element_blank()) +
                theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour="black"))
            plots[[i]] <- plot # add each plot into plot list
        }
    }

    return (plots)
}



#' plotPTMShiftOverallData
#'
#' This fucntion plots the shift on isoelectric point estimation produt of any modification (PTM) of all the predicted varaibles vs the expted variable (first column).
#' Both data1 as data2 must have the same size and shape.
#'
#' @param data1 the data.frame with estimations excluding PTM contribution
#' @param data2 the data.frame with estimations including PTM contribution
#'
#'

plotPTMShiftOverallData <- function (dat1, dat2){
    nm <- names(dat1)
    plots <- list()  # new empty list
    for (i in nm)
    {
        if(dat1[1L] != dat1[i])
        {
            #getting subset to anlysis. (x: expected variable, y1: predicted excluding PTM, y2: predicted including PTM)
            newData <- data.frame(x = dat1[1L], y1 = dat1[i], y2 = dat2[i])
            colnames(newData) <- c("x", "y1", "y2")

            #getting coefficients from lineal model.
            r_sqrt1 <- getCoeff(newData,newData$x, newData$y1)
            r_sqrt2 <- getCoeff(newData,newData$x, newData$y2)

            #reshape data to plot into groups (i.e. modified and non-modified).
            dataStats <- reshape(newData,
                                 varying= 2:3,
                                 v.names = "y",
                                 times=c("excluding modification","including modification"),
                                 direction = "long",
                                 sep = "")

            #print(fractions_stats)
            plot <- ggplot(dataStats, aes(x=x, y=y, shape=time))+
                geom_point(size=2.5, alpha=.4)+                          #set points size and transparency
                scale_shape_manual(values=c(1,2)) +                      #set points shape
                geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +       #add lineal model
                annotate("text", label=parseCoeff(r_sqrt1), parse=TRUE, x=Inf, y=-Inf, hjust=3.5, vjust= -5.0)+
                annotate("text", label=parseCoeff(r_sqrt2), parse=TRUE, x=Inf, y=-Inf, hjust=1.1, vjust= -.5)+
                #coord_fixed(ratio=1/2) +
                scale_y_continuous(breaks=seq(4, 10, 1)) +
                scale_x_continuous(breaks=seq(4, 7, .25))+
                xlab (nm[1L]) +
                ylab(i) +
                theme_bw()+
                theme(legend.position="none")+
                theme(legend.background=element_blank()) +
                theme(legend.key=element_blank()) +
                theme(panel.border = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour="black"))
            plots[[i]] <- plot # add each plot into plot list
        }
    }

    return (plots)
}



#' getCoeff
#'
#' This function retrive coefficient from a lineal model
#'
#' @param  data the data.frame containing data to anlysis
#' @param  xData the x variable to anlysis
#' @param  yData the y variable to anlysis
#'

getCoeff <- function(data, xData, yData){
    model <- lm(xData~yData, data)
    r2 = format(summary(model)$r.squared, digits=2)
    return (r2) #return correlation
}


#' parseCoeff
#'
#' This function parse numerical value into formule or expression
#'
#' @param  value the value to include into any expression
#'

parseCoeff <- function(coeff){
    eqn <- as.character(as.expression(
        substitute(italic(r)^2 == r2,
                   list(r2 = coeff))))
    return (eqn)
}
