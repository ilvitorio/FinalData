#Function Definitions and library requirements

#Libraries required

require("BMA")
require("PerformanceAnalytics")
require("ggplot2")
require("quantmod")
require("ROCR")

#User defined Functions

initialTransform <- function(data){
        #The function makes the initial data transformation from character to numeric
        numVar<-dim(data)
        data[,1]<-as.Date(data[,1],format="%d/%m/%Y")
        data[,-1]<-apply(data[,-1],2,as.numeric)
        return(data)
}

percentileRank <- function(vectorData){
        #The function ranks a generic ranking
        aux<-round(rank(vectorData)/length(vectorData),1)
        return(aux)
}

DFRank <- function(DFxts){
        #The function ranks the data frame and returns a DF of Ranking
        burn_date<-'1999-12'
        initial_index<-paste("/",burn_date,sep="")
        indexer<-paste(burn_date,"/",sep="")
        index_dates<- index(DFxts[indexer])[-1]
        DFRanked<-DFxts
        DFRanked[initial_index,] <- apply(DFxts[initial_index],2,percentileRank)
        
        for (i in index_dates) { 
                #Computes the roll over in the XTS object for taking into account the 
                #date up to the time
                newIndex<-paste("/",as.Date(i),sep="")
                DFtemp<-DFxts[newIndex]
                newDecile<-tail(apply(DFtemp,2,percentileRank),n=1)
                DFRanked[as.Date(i),]<-newDecile
        }
        
        return(DFRanked)
}


trendExtract <- function(tsUni,frequency){
        #Perform a STL decomposition and extrats the trending time series
        STLObj<-stl(ts(tsUni,frequency=frequency),s.window="periodic")
        TObj<-STLObj$time.series[,2]
        return(TObj)
        
}

DFTrend <- function(DFxts,frequency){
        #The function performs a Seasonal-Trend Decomposition on each column 
        #and returns the DF of the trending times series in a ts object
        DFTrendObj <- apply(DFxts,2,trendExtract,frequency=frequency)
        startDate<-as.Date(start(DFxts))
        Year<-as.numeric(format(startDate,"%Y"))
        Month<-as.numeric(format(startDate,"%m"))
        DFTrendObj <- ts(DFTrendObj,start=c(Year,Month),frequency=frequency)
        return(DFTrendObj)
}


tsMMCycle<- function(tsUni,radius){
        #This function returns a dataframe of the coordinates of the maximum and minimum 
        #points in a time series.
        
        #Take all the critical points in the time series
        diffSeries<- diff(sign(diff(tsUni)))
        points_max<-which(diffSeries==-2)+1
        points_min<-which(diffSeries==+2)+1
        
        #Filter the cycle meaningful critical points in the Time Series following a modified
        #Newton search approach
        auxlagged<- cbind(tsUni,
                          lag(tsUni,radius),
                          lag(tsUni,-radius))
        decision<- sign((auxlagged[,2]-auxlagged[,1])*(auxlagged[,1]-auxlagged[,3]))
        #Filtering the critical points into the Meaningful ones.
        points_max_clear<- points_max[decision[points_max+radius] == -1]
        points_min_clear<- points_min[decision[points_min+radius] == -1]
        
        #Construct the Iterator for filling the vector of expansion and contraction. 
        tsCycle <- rep(0,length(tsUni))
        first_point<-min(points_max,points_min)
        
        #Initialization of the Vector
        if (tsUni[1]<tsUni[first_point]) {
                tsCycle[1]<-0
        }else{
                tsCycle[1]<-1
        }
        
        #Fill the time series of 1s and 0s depending on the cycle.
        tsCycle[points_max_clear] <- rep(1,length(points_max_clear))
        tsCycle[points_min_clear] <- rep(-1,length(points_min_clear))
        tsCycle<- cumsum(tsCycle)
        
        return(tsCycle)
}

DFCycles <- function(DFTrend,radius){
        #THe function takes a smoothed time series dataframe and delivers a 1 and 0 for
        #increasing and decreasing periods.
        # 0 For increase
        # 1 For decrease
        # This is consistent with the NBER variable definition
        colnumber<-dim(DFTrend)[2]
        rownumber<-dim(DFTrend)[1]
        DFCycles<-DFTrend
        
        for (i in seq(colnumber) ){
                result<-try(tsMMCycle(DFTrend[,i],radius),silent=TRUE)
                if (class(result)=="try-error") {
                        DFCycles[,i]<-rep(NA,rownumber)
                } else{
                        DFCycles[,i]<-result
                }
                
        }
        
        startDate<-start(DFTrend)
        frequency<-tsp(DFTrend)[3]
        DFCycles <- ts(DFCycles,start=startDate,frequency=frequency)
        
        return(DFCycles)
        
}

DFStage <- function(tsDec,tsCycle,minT,maxT){
        #The function puts in place the Economic Bussines Cycles for a variable based
        #in its deciles and its cycle state.
        
        #Heritance of the time series variable and Default State
        tsStage<-tsDec
        tsStage<-rep("Equilibrium", length(tsStage))
        
        ### Filtering the Equilibrium to Peak
        index_threshold <- tsDec >= maxT  
        index_composed <- index_threshold  & !(as.logical(as.numeric(tsCycle)))
        tsStage[index_composed]<- c("Equilibrium to Peak")
        
        ### Filtering Peak to Equilibirum
        index_threshold <-  tsDec >= maxT
        index_composed <- index_threshold  & (as.logical(as.numeric(tsCycle)))
        tsStage[index_composed]<- c("Peak to Equilibrium")
        
        ### Filtering Trough to Equilibrium
        index_threshold <-  tsDec <= minT
        index_composed <- index_threshold  & !(as.logical(as.numeric(tsCycle)))
        tsStage[index_composed]<- c("Pit to Equilibrium")
        
        ### Filtering Equilibirum to Trough 
        index_threshold <-  tsDec <= minT
        index_composed <- index_threshold  & (as.logical(as.numeric(tsCycle)))
        tsStage[index_composed]<- c("Equilibrium to Pit")
        ##
        
        return(tsStage)
}

extractDate <- function(tsUni){
        #This function extract the date from a Time Series and leaves it
        #as a date format
        dates<-as.Date(index(as.xts(tsUni)))
        return(dates)
}

ts2xts <- function(tsUni){
        #Trasnform the ts object to an xts
        date<-extractDate(tsUni)
        xtsUni<-xts(tsUni, order.by =  date_test)
        return(xtsUni)
}

stringMonth2day <- function(xtsUniM,xtsDF){
        #Fill a Month Label to a daily timeseries.
        
        xtsADF <- xtsDF
        xtsADF$Cycle <- rep(NA)
        for (i in unique(xtsUniM)){
                string<-as.character(i)
                month_dates <- format(as.Date(index(xtsUniM[xtsUniM==string])),"%Y-%m")
                xtsADF[month_dates,"Cycle"]   <- c(string)
        }
        #xtsADF[,-dim(xtsADF)] <- as.numeric( xtsADF[,-dim(xtsADF)] )
        return(xtsADF)
}


DFsummarymean <- function(DFlabeled){
        #The function makes a summary of all the columns by only one label column
        col_names  <-names(DFlabeled)
        col_nums   <-length(col_names)
        order_index <-order(unique(DFlabeled[,col_nums]))
        row_names   <- unique(DFlabeled[,col_nums])[order_index]
        colDF_names   <- col_names[-c(1,col_nums)]
        DFsummary<-data.frame( matrix(NA,length(row_names),length(colDF_names) ))
        colnames(DFsummary) <- colDF_names
        rownames(DFsummary) <- row_names
        
        for (i in colDF_names ) {
                DFsummary[,i]<-aggregate(x = as.numeric(as.character(DFlabeled[,i])) , 
                                         by = as.list(DFlabeled[col_names[col_nums]]), 
                                         FUN = "mean")[,2]
        }
        
        
        return(DFsummary)
}


#DFmonthday <- function()

#Begins the coding

## Charge the data and make the initial data transformation
#Data_Month_Level_R   <- initialTransform(read.csv("~/General Documents/FED Study/Betalvl2-master/Datasets/Data_Month_Level_R.csv",stringsAsFactors=FALSE))
#Data_Month_Perc_R    <- initialTransform(read.csv("~/General Documents/FED Study/Betalvl2-master/Datasets/Data_Month_Perc_R.csv",stringsAsFactors=FALSE))
#Data_Month_PercCh_R  <- initialTransform(read.csv("~/General Documents/FED Study/Betalvl2-master/Datasets/Data_Month_PercCh_R.csv",stringsAsFactors=FALSE))
#Data_Quarter_Level_R <- initialTransform(read.csv("~/General Documents/FED Study/Betalvl2-master/Datasets/Data_Month_PercCh_R.csv",stringsAsFactors=FALSE))
#Data_Index_Level_R   <- initialTransform(read.csv("~/General Documents/FED Study/Betalvl2-master/Datasets/IndexDataR_Filter.csv",stringsAsFactors=FALSE))


##Chargin the data Second Source
Data_Month_Level_R   <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Data_Month_Level_R.csv",stringsAsFactors=FALSE))
Data_Month_Perc_R    <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Data_Month_Perc_R.csv",stringsAsFactors=FALSE))
Data_Month_PercCh_R  <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Data_Month_PercCh_R.csv",stringsAsFactors=FALSE))
Data_Quarter_Level_R <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/Data_Month_PercCh_R.csv",stringsAsFactors=FALSE))
Data_Index_Level_R   <- initialTransform(read.csv("~/Upwork/Betasmartz - Black Litterman/FED's Model/IndexDataR_Filter.csv",stringsAsFactors=FALSE))


## Transform the variables into XTS DataFrames
Data_Month_Level_ts   <- xts(Data_Month_Level_R[,-1] , order.by=Data_Month_Level_R[,1])
Data_Month_Perc_ts    <- xts(Data_Month_Perc_R[,-1] , order.by=Data_Month_Perc_R[,1],frequency=12)
Data_Month_PercCh_ts  <- xts(Data_Month_PercCh_R[,-1] , order.by=Data_Month_PercCh_R[,1])
Data_Quarter_Level_ts <- xts(Data_Quarter_Level_R[,-1] , order.by=Data_Quarter_Level_R[,1])
Data_Index_Level_ts   <- xts(Data_Index_Level_R[,-1] , order.by=Data_Index_Level_R[,1])

## Filter all the cases to work with
Data_Month_Level_clear   <- Data_Month_Level_ts[complete.cases(Data_Month_Level_ts)]
Data_Month_Perc_clear    <- Data_Month_Perc_ts[complete.cases(Data_Month_Perc_ts)]
Data_Month_PercCh_clear  <- Data_Month_PercCh_ts[complete.cases(Data_Month_PercCh_ts)]
Data_Quarter_Level_clear <- Data_Quarter_Level_ts[complete.cases(Data_Quarter_Level_ts)]
Data_Index_Level_Clear   <- Data_Index_Level_ts[complete.cases(Data_Index_Level_ts)]

## Apply logarithms
Data_Month_Level_log   <- log(Data_Month_Level_clear)
Data_Index_Level_log <- log(Data_Index_Level_Clear)

#Data_Quarter_Level_log <- log(Data_Quarter_Level_clear)

## Returns of the data 1-Month, 3-Month ,6-Month and 12-Month from Log
Data_Month_Level_R1 <- diff(Data_Month_Level_log,1)
Data_Month_Level_R3 <- diff(Data_Month_Level_log,3)
Data_Month_Level_R6 <- diff(Data_Month_Level_log,6)
Data_Month_Level_R12 <- diff(Data_Month_Level_log,12)

## Returns of the data 1-Month, 3-Month ,6-Month and 12-Month from Perc
Data_Month_Perc_R1 <- diff(Data_Month_Perc_clear,1)
Data_Month_Perc_R3 <- diff(Data_Month_Perc_clear,3) #Taking this one for FED Factor
Data_Month_Perc_R6 <- diff(Data_Month_Perc_clear,6)
Data_Month_Perc_R12 <- diff(Data_Month_Perc_clear,12)

## Returns of the data 1-Month, 3-Month, 6-Month and 12-Month from Index
Data_Index_Level_R1 <- diff(Data_Index_Level_log,1)
Data_Index_Level_R3 <- diff(Data_Index_Level_log,3)
Data_Index_Level_R6 <- diff(Data_Index_Level_log,6)
Data_Index_Level_R12 <- diff(Data_Index_Level_log,12)


## Filter again the cases to work with
Data_Month_Level_R1  <-Data_Month_Level_R1[complete.cases(Data_Month_Level_R1)]
Data_Month_Level_R3  <-Data_Month_Level_R3[complete.cases(Data_Month_Level_R3)]
Data_Month_Level_R6  <-Data_Month_Level_R6[complete.cases(Data_Month_Level_R6)]
Data_Month_Level_R12 <-Data_Month_Level_R12[complete.cases(Data_Month_Level_R12)]


Data_Month_Perc_R1  <-Data_Month_Perc_R1[complete.cases(Data_Month_Perc_R1)]
Data_Month_Perc_R3  <-Data_Month_Perc_R3[complete.cases(Data_Month_Perc_R3)]
Data_Month_Perc_R6  <-Data_Month_Perc_R6[complete.cases(Data_Month_Perc_R6)]
Data_Month_Perc_R12 <-Data_Month_Perc_R12[complete.cases(Data_Month_Perc_R12)]

Data_Index_Level_R1  <-Data_Index_Level_R1[complete.cases(Data_Index_Level_R1)]
Data_Index_Level_R3  <-Data_Index_Level_R3[complete.cases(Data_Index_Level_R3)]
Data_Index_Level_R6  <-Data_Index_Level_R6[complete.cases(Data_Index_Level_R6)]
Data_Index_Level_R12 <-Data_Index_Level_R12[complete.cases(Data_Index_Level_R12)]


## Rank the data by deciles

#Returns by Deciles
Data_Month_Level_R1_Dec <- DFRank(Data_Month_Level_R1)
Data_Month_Level_R3_Dec <- DFRank(Data_Month_Level_R3)
Data_Month_Level_R6_Dec <- DFRank(Data_Month_Level_R6)
Data_Month_Level_R12_Dec<- DFRank(Data_Month_Level_R12)

Data_Month_Perc_R1_Dec <- DFRank(Data_Month_Perc_R1)
Data_Month_Perc_R3_Dec <- DFRank(Data_Month_Perc_R3)
Data_Month_Perc_R6_Dec <- DFRank(Data_Month_Perc_R6)
Data_Month_Perc_R12_Dec <- DFRank(Data_Month_Perc_R12)

#Percentage Variables by Deciles
Data_Month_Perc_Dec    <- DFRank(Data_Month_Perc_clear)
Data_Month_PercCh_Dec  <- DFRank(Data_Month_PercCh_clear)

#Seasonal Trend Decomposition to find local minima and maxima
Data_Month_Perc_Trend <- DFTrend(Data_Month_Perc_Dec,12)

#Cycle DataFrame
Data_Month_Perc_Cycle <- DFCycles(Data_Month_Perc_Trend,12)


#Stage DataFrame
minimumT<-0.5
maximumT<-0.6
CycleVar <- DFStage(Data_Month_Perc_Trend[,24],Data_Month_Perc_Cycle[,24],minimumT,maximumT)

dataplot<-data.frame(cbind(Data_Month_Perc_Trend[,24]),CycleVar)
colnames(dataplot)<-c("a","b")
dataplot$a<-as.numeric(as.character(dataplot$a))
line<-ggplot(dataplot, aes(x=seq(dim(dataplot)[1]), y=a, color=b)) + ggtitle("Economic Cycle") + xlab("Time") + ylab("Output GDP GAP") + geom_point()

#Extract the date and coerce it 
date_test<-extractDate(Data_Month_Perc_Trend[,24]) -1
CycleVar_ts <- xts(CycleVar, order.by =  date_test)
colnames(CycleVar_ts) <- c("Cycle")

#Create the Index - Bussiness Cycle DataFrame
Index_R1_Cycle_ts  <- stringMonth2day(CycleVar_ts,Data_Index_Level_R1)
#Filtering the non available data
Index_R1_Cycle_ts  <- Index_R1_Cycle_ts[complete.cases(Index_R1_Cycle_ts)]

#Transform the xts matrix into a dataframe for input extraction
Index_Cycle        <- data.frame(date=index(Index_R1_Cycle_ts) , coredata(Index_R1_Cycle_ts))

#Expected Return Summary Data by Stage 
testDFMean<- DFsummarymean(Index_Cycle)


### Building the BMA Model 

#The Data comes from here
Data_Month_Level_R3
Data_Month_Perc_R3

# Data_Month_PercCh_clear not inclided in this, due to their short length
dummyCycle <- xts(model.matrix( ~Cycle-1,data=CycleVar_ts), order.by = index(CycleVar_ts))
lag12_dummyCycle <- lag(dummyCycle,-12)
lag12_dummyCycle <- lag12_dummyCycle[complete.cases(lag12_dummyCycle)]

#Total data DataFrame 
Total_Data <- merge(Data_Month_Level_R3,Data_Month_Perc_R3,dummyCycle)
#CLean the total DataFrame
Total_Data<-Total_Data[complete.cases(Total_Data)]

lag12_Data <- merge(Data_Month_Level_R3,Data_Month_Perc_R3,lag12_dummyCycle)
#Clean the xts object
lag12_Data<-lag12_Data[complete.cases(lag12_Data)]

#XTS to DF
lag12_Data_DF <- data.frame(lag12_Data)


### Beggining the BMA Adjusting


## Equilibrium adjust
formula_Equi <-as.formula(paste(colnames(lag12_Data)[54], 
                               "~",
                               paste(colnames(lag12_Data)[seq(43)], 
                                     collapse = "+"),
                                sep = ""))
        

## Adjust the model and Store it
BMA.Equilibrium <- bic.glm(formula_Equi,data=lag12_Data_DF, glm.family=binomial(link="probit"), OR=500, OR.fix=50,nBest=250 )
imageplot.bma(BMA.Equilibrium,order="probne0")


## Equilibrium to Peak adjust
formula_Equi2Peak <-as.formula(paste(colnames(lag12_Data)[55], 
                                "~",
                                paste(colnames(lag12_Data)[seq(43)], 
                                      collapse = "+"),
                                sep = ""))


## Adjust the model and Store it
BMA.Equilibrium2Peak <- bic.glm(formula_Equi2Peak,data=lag12_Data_DF, glm.family=binomial(link="probit"), OR=500, OR.fix=50,nBest=250 )
imageplot.bma(BMA.Equilibrium2Peak,order="probne0")


## Equilibrium to Pit adjust
formula_Equi2Pit <-as.formula(paste(colnames(lag12_Data)[56], 
                                     "~",
                                     paste(colnames(lag12_Data)[seq(43)], 
                                           collapse = "+"),
                                     sep = ""))


## Adjust the model and Store it
BMA.Equilibrium2Pit <- bic.glm(formula_Equi2Pit,data=lag12_Data_DF, glm.family=binomial(link="probit"), OR=500, OR.fix=50,nBest=250 )
imageplot.bma(BMA.Equilibrium2Pit,order="probne0")



## Peak to Equilibrium adjust
formula_Peak2Equi <-as.formula(paste(colnames(lag12_Data)[57], 
                                    "~",
                                    paste(colnames(lag12_Data)[seq(43)], 
                                          collapse = "+"),
                                    sep = ""))


## Adjust the model and Store it
BMA.Peak2Equilibrium <- bic.glm(formula_Peak2Equi,data=lag12_Data_DF, glm.family=binomial(link="probit"), OR=500, OR.fix=50,nBest=250 )
imageplot.bma(BMA.Peak2Equilibrium,order="probne0")


## Pit to Equilibrium adjust
formula_Pit2Equi <-as.formula(paste(colnames(lag12_Data)[58],
                                     "~",
                                     paste(colnames(lag12_Data)[c(seq(17),20:43)], 
                                           collapse = "+"),
                                     sep = ""))


## Adjust the model and Store it
BMA.Pit2Equilibrium <- bic.glm(formula_Pit2Equi,data=lag12_Data_DF, glm.family=binomial(link="probit"), OR=500, OR.fix=50,nBest=250 )
imageplot.bma(BMA.Pit2Equilibrium,order="probne0")


### Proability predictions for 12 month ahead.

predict_Equilibrium <- predict(BMA.Equilibrium, newdata = Total_Data)
predict_Equilibrium2Peak <- predict(BMA.Equilibrium2Peak, newdata = Total_Data)
predict_Equilibrium2Pit <- predict(BMA.Equilibrium2Pit, newdata = Total_Data)
predict_Peak2Equilibrium <- predict(BMA.Peak2Equilibrium, newdata = Total_Data)
predict_Pit2Equilibrium <- predict(BMA.Pit2Equilibrium, newdata = Total_Data)


DF_12M_P <- cbind(predict_Equilibrium,predict_Equilibrium2Peak,predict_Equilibrium2Pit,predict_Peak2Equilibrium,predict_Pit2Equilibrium)


#Covariance matrices

#Equilibrium Covariance
Equi_Matrix<-Index_Cycle[Index_Cycle["Cycle"]=="Equilibrium",-c(1,13)]
Equi_Matrix[, seq(11)] <- sapply(Equi_Matrix[,seq(11)], as.character)
Equi_Matrix[, seq(11)] <- sapply(Equi_Matrix[,seq(11)], as.numeric)
Cov_Equi <- cov(Equi_Matrix)

#Equilibrium to Peak Covariance
Equi2Peak_Matrix<-Index_Cycle[Index_Cycle["Cycle"]=="Equilibrium to Peak",-c(1,13)]
Equi2Peak_Matrix[, seq(11)] <- sapply(Equi2Peak_Matrix[,seq(11)], as.character)
Equi2Peak_Matrix[, seq(11)] <- sapply(Equi2Peak_Matrix[,seq(11)], as.numeric)
Cov_Equi2Peak <- cov(Equi2Peak_Matrix)

#Pit to Equilibrium Covariance
Pit2Equi_Matrix<-Index_Cycle[Index_Cycle["Cycle"]=="Pit to Equilibrium",-c(1,13)]
Pit2Equi_Matrix[, seq(11)] <- sapply(Pit2Equi_Matrix[,seq(11)], as.character)
Pit2Equi_Matrix[, seq(11)] <- sapply(Pit2Equi_Matrix[,seq(11)], as.numeric)
Cov_Pit2Equi <- cov(Pit2Equi_Matrix)

#Equilibrium to Pit Covariance
Equi2Pit_Matrix<-Index_Cycle[Index_Cycle["Cycle"]=="Equilibrium to Pit",-c(1,13)]
Equi2Pit_Matrix[, seq(11)] <- sapply(Equi2Pit_Matrix[,seq(11)], as.character)
Equi2Pit_Matrix[, seq(11)] <- sapply(Equi2Pit_Matrix[,seq(11)], as.numeric)
Cov_Equi2Pit <- cov(Equi2Pit_Matrix)


# Output Excel writing

#Covariance Matrices
write.csv(Cov_Equi2Pit, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/Cov_Equi2Pit.csv")
write.csv(Cov_Pit2Equi, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/Cov_Pit2Equi.csv")
write.csv(Cov_Equi, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/Cov_Equi.csv")
write.csv(Cov_Peak2Equi, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/Cov_Peak2Equi.csv")
write.csv(Cov_Equi2Peak, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/Cov_Equi2Peak.csv")

#Expected Returns
write.csv(testDFMean, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/testDFMean.csv")

#Probabilities for each model
write.csv(DF_12M_P, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/predict_Probs12.csv")

#Index Cycle Variable
write.csv(Index_Cycle, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/Index_Cycle.csv")

#Write the Total Data
write.csv(Total_Data, "C:/Users/Vitty2/Documents/Upwork/Betasmartz - Black Litterman/FED's Model/Total_Data.csv")

# a<-ts(as.vector(Data_Month_Perc_Dec$ROUTGAP),start=as.Date(start(Data_Month_Perc_Dec$ROUTGAP)),frequency=12)
# stl(a,s.window="periodic")
# 
# radius<-12
# lagged<-cbind(Data_Month_Perc_Trend[,24],lag(Data_Month_Perc_Trend[,24],radius),lag(Data_Month_Perc_Trend[,24],-radius))
# colnames(lagged)<-c("orig","atras","adela")
# testiviri<-cbind(lagged[,2]-lagged[,1],lagged[,1]-lagged[,3],lagged[,1],diff(sign(diff(Data_Month_Perc_Trend[,24]))))
# 
# colnames(testiviri)<- c("var1","var2","orig","detector")