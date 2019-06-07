setGeneric(name = "chromatogramPlot", 
           def = function(object,
                          title.size = 15,
                          lab.size = 15,
                          axis.text.size =15,
                          alpha = 0.5,
                          title = ""){
            info <- object@phenoData@data
            data <- object@.Data
            rm(list = c("object"))
            data <- apply(data, 2, function(x){
              x <- x[[1]]
              x <- data.frame("mz" = x@rtime, "intensity" = x@intensity, 
                              stringsAsFactors = FALSE)
              list(x)
            })
            
            data <- lapply(data, function(x){
              x[[1]]
            })
            
            data <- mapply(FUN = function(x, y){
              x <- data.frame(x, "group" = y, stringsAsFactors = FALSE)
              list(x)
            },
            x = data,
            y = info[,2])
            
            data <- do.call(rbind, args = data)
           
            plot <- ggplot2::ggplot(data = data, 
                                    ggplot2::aes(x = mz, y = intensity)) +
              ggplot2::geom_line(data = data, mapping = ggplot2::aes(colour = group), alpha = alpha) +
              ggplot2::theme_bw() +
              ggplot2::labs(x = "Mass to charge ratio (m/z)", y = "Intensity", title = title) +
              ggplot2::theme(
                plot.title = ggplot2::element_text(color="black", size = title.size,
                                                   face = "plain",
                                                   hjust = 0.5),
                axis.title = ggplot2::element_text(color="black", size = lab.size,
                                                   face = "plain"),
                axis.text = ggplot2::element_text(color="black", size = axis.text.size,
                                                  face = "plain")
              )
             
            
           })
 


setGeneric(name = "plotAdjustedRT",
           def = function(object,
                          title.size = 15,
                          lab.size = 15,
                          axis.text.size =15){
             diffRt <- xcms::rtime(object, adjusted = TRUE) - xcms::rtime(object,
                                                              adjusted = FALSE)
             diffRt <- split(diffRt, MSnbase::fromFile(object))
             xRt <- xcms::rtime(object,
                                adjusted = TRUE,
                                bySample = TRUE)
             
             sample_name <- object@phenoData@data$sample_name
             sample_group <- object@phenoData@data$sample_group
             
             diffRt <- mapply(FUN = function(x, y){
               list(data.frame(x, y, stringsAsFactors = FALSE))
             },
             x = diffRt,
             y = sample_name)
           
             xRt <- mapply(FUN = function(x, y){
               list(data.frame(x, y, stringsAsFactors = FALSE))
             },
             x = xRt,
             y = sample_name)
           
             diffRt <- do.call(what = rbind, args = diffRt)
            xRt <- do.call(rbind, xRt)     
           
            temp.data <- data.frame(xRt, diffRt, stringsAsFactors = FALSE)
           
            colnames(temp.data) <- c("rt", "sample.name","diffRT", "sample.name2") 
           rm(list = c("object", "xRt", "diffRt")) 
           
           plot <- ggplot2::ggplot(data = temp.data, ggplot2::aes(x = rt, y = diffRT)) +
             ggplot2::geom_line(data = temp.data, ggplot2::aes(color = sample.name)) + 
             ggplot2::theme_bw() + 
             ggplot2::labs(x = "Retention time (second)", y = "RT deviation (second)") +
             ggplot2::theme(legend.position = "none",
                            axis.title = ggplot2::element_text(color="black", size = lab.size,
                                                               face = "plain"),
                            axis.text = ggplot2::element_text(color="black", size = axis.text.size,
                                                              face = "plain"))
           
           plot
           
              
           })













