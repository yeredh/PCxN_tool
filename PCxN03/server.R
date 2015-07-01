
ExtractTopConnected <- function(path.dframe,path.target,top.n=25,topMethod="abs",p.cut=0.05,cor.cut=0.25){
  # Function to return the top n pathways connected with the target pathway
  # 
  # Args:
  #   path.dframe: a data frame with the names of the pathways
  #                (Pathway.A and Pathway.B) and their corresponding
  #                value for the pathway correlation (PathCor)
  #   path.target: string with the name of the pathway of interest
  #   top.n: the cut-off for the top n pathways connected to the 
  #          the pathway of interest.
  #   topMethod: method to sort the edge weights (PathCor). There are
  #              3 methods available
  #              1. "abs" sort by the absolute value of the edge weight
  #                 (|PathCor|)
  #              2. "decreasing" sort the edge weight in decreasing order
  #              2. "increasing" sort the edge weight in increasing order
  #
  # Returns:
  #   A data frame with the same format as the input data frame path.dframe
  #   with the top n pathways connected to the pathway of interest
  
  # applu cut-offs
  path.dframe <- path.dframe[p.Adjust < p.cut & abs(PathCor) > cor.cut,]
  
  # get all edges for target pathway
  setkey(path.dframe,Pathway.A)
  path.subnet = path.dframe[path.target,]
  setkey(path.dframe,Pathway.B)
  path.subnet = rbind(path.subnet,path.dframe[path.target,])
  
  # get top edges based on type of ordering
  selectMethod <- c("abs","decreasing","increasing") %in% topMethod
  if(selectMethod[1]){
    # sort by |PathCor|
    path.subnet <- path.subnet[order(abs(PathCor),decreasing=T),]
  }
  if(selectMethod[2]){
    path.subnet <- path.subnet[order(PathCor,decreasing=T),]
  }
  if(selectMethod[3]){
    path.subnet <- path.subnet[order(PathCor,decreasing=F),]
  }
  if(sum(selectMethod)==0){
    print("Not a valid choice to sort edges")
    return()
  }
  # return top edges
  if(dim(path.subnet)[1] >= top.n){
    return(path.subnet[1:top.n,])
  }else{
    return(path.subnet)
  }
  
}

SwapPathwayNames <- function(path.dframe,path.target,top.n){
  # Auxiliary function to swap the pathways names in a data frame so that
  # the target pathway name is in the Pathway.A and the top n connected
  # pathway names are in the Pathway.B column
  #
  # Args:
  #   path.dframe: a data frame with the names of the pathways
  #                (Pathway.A and Pathway.B) 
  #   path.target: string with the name of the pathway of interest
  #   top.n: the cut-off for the top n pathways connected to the 
  #          the pathway of interest.
  #
  # Returns:
  #   A data frame with the same format as the input data frame path.dframe
  
  swap.ind <- which(path.dframe$Pathway.B[1:top.n] == path.target)
  path.dframe$Pathway.B[swap.ind] <- path.dframe$Pathway.A[swap.ind]
  path.dframe$Pathway.A[swap.ind] <- path.target
  return(path.dframe)
}


ExtractSubNet2 <- function(path.dframe,path.target,top.n=25,top.k=NULL,topMethod="abs",
                           p.cut=0.05,cor.cut=0.25){
  # Function to return the top n pathways connected with the target pathway and the
  # edges among the top n pathways connected with the pathway of interest
  #
  # Args:
  #   path.dframe: a data frame with the names of the pathways
  #                (Pathway.A and Pathway.B) and their corresponding
  #                value for the pathway correlation (PathCor)
  #   path.target: string with the name of the pathway of interest
  #   top.n: the cut-off for the top n pathways connected to the 
  #          the pathway of interest.
  #   top.k: the cut-off for the number of edges to include among the
  #          top connected pathways. If top.n=NULL, then all connections
  #          between the top pathways are considered
  #   topMethod: method to sort the edge weights (PathCor). There are
  #              3 methods available
  #              1. "abs" sort by the absolute value of the edge weight
  #                 (|PathCor|)
  #              2. "decreasing" sort the edge weight in decreasing order
  #              2. "increasing" sort the edge weight in increasing order
  #
  # Returns:
  #   A data frame with the same format as the input data frame path.dframe
  #   with the top n pathways connected to the pathway of interest and the 
  #   connections among the top connected pathways
  
  
  subnet.dframe <- ExtractTopConnected(path.dframe,path.target,top.n,topMethod,p.cut,cor.cut)
  if(dim(subnet.dframe)[1] == 0){
    cat("No pathways found! >=[ \n")
    fup = rep(NA,dim(subnet.dframe)[2])
    names(fup) = colnames(subnet.dframe)
    return(as.data.frame(t(fup)))
  }
  if(dim(subnet.dframe)[1] < top.n){
    top.n <- dim(subnet.dframe)[1]
  }
  # swap pathway names
  subnet.dframe <- SwapPathwayNames(subnet.dframe,path.target,top.n)
  return(subnet.dframe)
}


TopTable2 <- function(path.dframe,top.n,path.target){
  if(all(is.na(path.dframe))){
    cat("No pathways found! >=[ \n")
    fup = rep(NA,dim(path.dframe)[2])
    names(fup) = colnames(path.dframe)
    return(as.data.frame(t(fup)))
  }
  path.ind <- c(grep(path.target,path.dframe$Pathway.A,fixed=T),
                grep(path.target,path.dframe$Pathway.B,fixed=T))
  if(length(path.ind) < top.n){
    top.n <- length(path.ind)
  }
  path.dframe <- path.dframe[1:top.n,]
  
  
  path.dframe <- SwapPathwayNames(path.dframe,path.target,top.n)
  tab <- data.frame(Pathway=rep(NA,top.n),
                    PathCor=rep(NA,top.n),
                    Overlap.Coefficient=rep(NA,top.n),
                    p.value=rep(NA,top.n))
  
  tab$Pathway <- path.dframe$Pathway.B
  tab$PathCor <- path.dframe$PathCor
  tab$Overlap.Coefficient <- path.dframe$Overlap.Coefficient
  tab$p.value <- path.dframe$p.Adjust
  return(tab)
}

library(shiny)

shinyServer(function(input, output){
  output$caption <- renderText({
    input$tool1
  })
  
  datasetInput <- reactive({
    sMethod <- switch(input$method,
                       "Absolute Value" = "abs",
                       "Decreasing" = "decreasing",
                       "Increasing" = "increasing")
    na.omit(ExtractSubNet2(path.dframe=path.net.dframe,
                  path.target=input$tool1,
                  top.n=input$top.n,
                  top.k=5,
                  topMethod=sMethod,
                  p.cut=input$p.cut,
                  cor.cut=input$cor.cut))
  })
  
  output$view <- renderTable({
    TopTable2(datasetInput(),input$top.n,input$tool1)
  },digits=4)
  

})