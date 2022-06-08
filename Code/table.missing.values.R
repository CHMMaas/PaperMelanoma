missings.table <- function(data){
  missing.values.table <- data.frame(m=rep(0, ncol(data)))
  rownames(missing.values.table) <- colnames(data)
  for (variable in colnames(data)){
    missing.values.table[variable,] <- sum(is.na(data[,variable]))
  }
  print(missing.values.table[order(missing.values.table$m, decreasing=TRUE), , drop=FALSE])  
}