code_poly <- function(df=NULL, predictor=NULL, poly.order=NULL, orthogonal=TRUE, draw.poly=TRUE){
  require(reshape2)
  require(ggplot2)
  # Codes raw or orthogonal polynomial transformations of a predictor variable
  # be sure to not have an actual variable named "predictor" in your data.frame
  # ultimately adds 2 or more columns, including:
  # (predictor).Index, and a column for each order of the polynomials 
  #
  # Written by Matt Winn: http://www.mattwinn.com/tools/R_add_polynomials_to_df.html
  # Revised by Dan Mirman (11/3/2014): 
  #   - call to poly now uses predictor.index and 1:max instead of unique() to deal with out-of-order time bins
  #   - combined indexing and alignment with original data, mainly to avoid problems when poly.order==1
  # Revised by Matt Winn (2/12/2015)
  #   - uses `[` instead of `$` to use variable name to extract column
  #   - computes polynomial on unique(sort(predictor.vector))
  #       rather than 1:max(predictor.vector) 
  #       to accomodate non-integer predictor (e.g. time) levels
  #   - Accomodates missing/unevenly-spaced time bins 
  #       by indexing each sorted unique time bin and using the index to extract
  #       the polynomial value
  
  #===========================================================#
  # convert choice for orthogonal into choice for raw
  raw <- (orthogonal-1)^2
  
  # make sure that the declared predictor is actually present in the data.frame
  if (!predictor %in% names(df)){
    warning(paste0(predictor, " is not a variable in your data frame. Check spelling and try again"))
  }
  
  # Extract the vector to be used as the predictor
  # predictor.vector <- df[,which(colnames(df)==predictor)]
  predictor.vector <- df[[which(colnames(df)==predictor)]]
  
  # create index of predictor (e.g. numbered time bins)
  # the index of the time bin will be used later as an index to call the time sample
  predictor.indices <- as.numeric(as.factor(predictor.vector))
  
  df$temp.predictor.index <- predictor.indices
  
  #create x-order order polys (orthogonal if not raw)
  predictor.polynomial <- stats::poly(x = unique(sort(predictor.vector)),  
                               degree = poly.order, raw=raw)
  
  # use predictor index as index to align 
  # polynomial-transformed predictor values with original dataset
  # (as many as called for by the polynomial order)
  df[, paste("poly", 1:poly.order, sep="")] <- 
    predictor.polynomial[predictor.indices, 1:poly.order]
  
  # draw a plot of the polynomial transformations, if desired
  if (draw.poly == TRUE){
    # extract the polynomials from the df
    df.poly <- unique(df[c(predictor, paste("poly", 1:poly.order, sep=""))])
    
    # melt from wide to long format
    df.poly.melt <- melt(df.poly, id.vars=predictor)
    
    # Make level names intuitive
    # don't bother with anything above 6th order. 
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly1"] <- "Linear"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly2"] <- "Quadratic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly3"] <- "Cubic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly4"] <- "Quartic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly5"] <- "Quintic"
    levels(df.poly.melt$variable)[levels(df.poly.melt$variable)=="poly6"] <- "Sextic"
    
    # change some column names for the output
    colnames(df.poly.melt)[colnames(df.poly.melt) == "variable"] <- "Order"
    
    poly.plot <- ggplot(df.poly.melt, aes(y=value, color=Order))+
      aes_string(x=predictor)+
      geom_line()+
      xlab(paste0(predictor, " (transformed polynomials)"))+
      ylab("Transformed value")+
      scale_color_brewer(palette="Set1")+
      theme_bw()
    
    print(poly.plot)
  }
  
  # restore correct column names
  colnames(df)[colnames(df) == "temp.predictor.index"] <- paste0(predictor,".Index")
  return(df)
}