
##------------------------------------------------------------------------------
#' @title getIdentificationTable
#' @description Get identification table from a metIdentifyClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param ... One or multiple metIdentifyClass objects.
#' @param candidate.num The number of candidates.
#' @param type The type of identification table. 
#' @return A identification table (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found 
#' https://jaspershen.github.io/metIdentify/articles/metIdentify.html
setGeneric(name = "getIdentificationTable", 
           def = function(...,
                          candidate.num = 3,
                          type = c("old", "new")
           ){
             candidate.num <- round(candidate.num)
             if(candidate.num <= 0) candidate.num <- 1
             
             object <- list(...)
             if (any(unique(unlist(lapply(object, class))) != "metIdentifyClass")) {
               stop("Only for metIdentifyClass\n")
             }
             
             type <- match.arg(type)
             ##if object number is 1
             if(length(object) == 1){
               object <- object[[1]]
               database <- object@database
               identification.result <- object@identification.result
               if(is.null(identification.result[[1]])){
                 return(NULL)
               }
               ##add database information
               identification.result <- lapply(identification.result, function(x){
                 if(nrow(x) > candidate.num){
                   x <- x[1:candidate.num,,drop = FALSE]
                 }
                 data.frame(x, "Database" = object@database, stringsAsFactors = FALSE)
               })
               
               peak.table <- object@ms1.data
               match.result <- object@match.result
               
               if(type == "old"){
                 identification.table <- as.data.frame(matrix(nrow = nrow(peak.table), ncol = 3))
                 colnames(identification.table) <- c("MS2.spectrum.name","Candidate.number", "Identification")
                 identification.table[match.result[,1],1] <- object@ms1.info$name[match.result[,2]]
                 item <- colnames(identification.result[[1]])
                 identification.result <- lapply(identification.result, function(x){
                   paste(apply(x, 1, function(y){
                     paste(paste(item, as.character(y), sep = ":"), collapse = ";")
                   }), collapse = "{}")
                 })
                 
                 identification.table$Identification[match(names(identification.result), identification.table$MS2.spectrum.name)] <- 
                   unlist(identification.result)
                 
                 identification.table$Candidate.number <- sapply(identification.table$Identification, function(x){
                   if(is.na(x)) return(0)
                   return(length(stringr::str_split(string = x, pattern = "\\{\\}")[[1]]))
                 })
                 identification.table <- data.frame(peak.table, identification.table, stringsAsFactors = FALSE)
               }else{
                 identification.table <- vector(mode = "list", length = nrow(peak.table))
                 names(identification.table)[match.result[,1]] <- object@ms1.info$name[match.result[,2]]
                 identification.table[match(names(identification.result), names(identification.table))] <- identification.result
                 peak.table <- apply(peak.table, 1, list)
                 peak.table <- lapply(peak.table, unlist)
                 
                 identification.table <- mapply(FUN = function(x, y){
                   if(all(is.na(y))){
                     temp <- as.data.frame(matrix(c(x, rep(NA, 14)), nrow = 1), stringsAsFactors = FALSE)
                     colnames(temp) <- c(names(x), c("Compound.name", "CAS.ID", "HMDB.ID",
                                                     "KEGG.ID", "Lab.ID", "Adduct", "mz.error", "mz.match.score",
                                                     "RT.error", "RT.match.score",
                                                     "CE", "SS", "Total.score", "Database"))
                     list(temp)
                   }else{
                     temp <- as.data.frame(matrix(rep(x, nrow(y)), nrow = nrow(y), byrow = TRUE), 
                                           stringsAsFactors = FALSE)
                     if(nrow(temp) > 1){
                       temp[2:nrow(temp), 2:ncol(temp)] <- ""
                     }
                     colnames(temp) <- names(x)
                     temp <- data.frame(temp, y, stringsAsFactors = FALSE)
                     list(temp)
                   }
                 },
                 x = peak.table,
                 y = identification.table)
                 
                 identification.table <- as.data.frame(do.call(rbind, identification.table))
               }
               rownames(identification.table) <- NULL
               return(identification.table)
             }else{
               ###if object number is not 1
               peak.table <- object[[1]]@ms1.data
               database.name <- unlist(lapply(object, function(x) x@database))
               
               identification.result <-
                 lapply(object, function(x) {
                   iden.result <- x@identification.result
                   iden.result <- mapply(function(x, y){
                     if(nrow(x) > candidate.num){
                       x <- x[1:candidate.num,,drop = FALSE]
                     }
                     x <- data.frame(x, "MS2.spectra.name" = y, 
                                     stringsAsFactors = FALSE)
                     list(x)
                   }, 
                   x = iden.result,
                   y = names(iden.result))
                   match.result <- x@match.result
                   # names(iden.result) <-
                   #   match.result$MS1.peak.name[match(names(iden.result), match.result$MS2.spectra.name)]
                   iden.result <-
                     do.call(rbind, iden.result)
                   Peak.name <- match.result$MS1.peak.name[match(iden.result$MS2.spectra.name, match.result$MS2.spectra.name)]
                   iden.result <-
                     data.frame(Peak.name, iden.result, stringsAsFactors = FALSE)
                   rownames(iden.result) <- NULL
                   iden.result
                 })
               
               identification.result <- mapply(function(x, y, z){
                 list(data.frame(x, "Database" = y, "Dataset.index" = z,
                                 stringsAsFactors = FALSE))
               },
               x = identification.result,
               y = database.name,
               z = 1:length(identification.result))
               
               identification.result <- do.call(rbind, identification.result)
               identification.result <- as.data.frame(identification.result)
               
               identification.result <- plyr::dlply(.data = identification.result, 
                                                    .variables = plyr::.(Peak.name), 
                                                    .fun = function(x){
                                                      x <- x[order(x$SS,decreasing = TRUE),,drop = FALSE]
                                                      if(nrow(x) > candidate.num){x <- x[1:candidate.num,,drop = FALSE]}
                                                      x
                                                    })
               
               # if(type == "old"){
               identification.table <- as.data.frame(matrix(nrow = nrow(peak.table), ncol = 3))
               colnames(identification.table) <- c("MS2.spectrum.name","Candidate.number", "Identification")
               # identification.table[match.result[,1],1] <- object@ms1.info$name[match.result[,2]]
               item <- colnames(identification.result[[1]])[-1]
               identification.result <- lapply(identification.result, function(x){
                 paste(apply(x[,-1], 1, function(y){
                   paste(paste(item, as.character(y), sep = ":"), collapse = ";")
                 }), collapse = "{}")
               })
               
               identification.table$Identification[match(names(identification.result), peak.table$name)] <- 
                 unlist(identification.result)
               
               identification.table$Candidate.number <- sapply(identification.table$Identification, function(x){
                 if(is.na(x)) return(0)
                 return(length(stringr::str_split(string = x, pattern = "\\{\\}")[[1]]))
               })
               
               identification.table$MS2.spectrum.name <- unlist(lapply(identification.table$Identification, function(x){
                 if(is.na(x)) return(NA)
                 temp.name <- stringr::str_extract_all(string = x, pattern = "MS2.spectra.name:mz[0-9\\.\\-\\_]{1,30}rt[0-9\\.\\-\\_]{1,30}")[[1]]
                 temp.name <- stringr::str_replace_all(string = temp.name, pattern = "MS2.spectra.name:", replacement = "")
                 paste(temp.name, collapse = ";")
               }))
               
               identification.table <- data.frame(peak.table, identification.table, stringsAsFactors = FALSE)
               # }else{
               #   
               # }
               return(identification.table)
             }
           })







##------------------------------------------------------------------------------
#' @title trans2newStyle
#' @description Transform old style identification table to new style..
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param identification.table Identification table from getIdentificationTable or getIdentificationTable2.
#' @return A identification table (data.frame).
#' @export
#' @seealso The example and demo data of this function can be found 
#' https://jaspershen.github.io/metIdentify/articles/metIdentify.html

setGeneric(name = "trans2newStyle", 
           def = function(identification.table){
             if(all(colnames(identification.table) != "Identification")){
               return(identification.table)
             }else{
               identification <- identification.table$Identification
               identification <- 
                 pbapply::pblapply(identification, function(x){
                   if(is.na(x)){
                     return(
                       data.frame("Compound.name" = NA, 
                                  "CAS.ID" = NA, 
                                  "HMDB.ID" = NA,
                                  "KEGG.ID" = NA,
                                  "Lab.ID" = NA,
                                  "Adduct" = NA,
                                  "mz.error" = NA,
                                  "mz.match.score" = NA,
                                  "RT.error" = NA,
                                  "RT.match.score" = NA,
                                  "CE" = NA,
                                  "SS" = NA,
                                  "Total.score" = NA,
                                  "Database" = NA, 
                                  stringsAsFactors = FALSE)
                     )
                   }else{
                     x <- stringr::str_split(string = x, pattern = "\\{\\}")[[1]][1]
                     Compound.name <- stringr::str_split(string = x, pattern = "\\;CAS\\.ID", n = 2)[[1]][1]
                     require(magrittr)
                     Compound.name <- Compound.name %>% 
                       stringr::str_replace(string = ., pattern = "Compound.name\\:", replacement = "")
                     
                     x <- stringr::str_split(string = x, pattern = ";")[[1]]
                     
                     x <- sapply(c("CAS.ID\\:", "HMDB.ID\\:", "KEGG.ID\\:", "Lab.ID\\:",
                                   "Adduct\\:", "mz.error\\:", "mz.match.score\\:",
                                   "RT.error\\:", "RT.match.score\\:",
                                   "CE\\:", "SS\\:", "Total.score\\:", "Database\\:"), function(y){
                                     y <- grep(y, x, value = TRUE) %>% 
                                       stringr::str_replace(string = ., pattern = paste(y,"\\:", sep = ""), 
                                                            replacement = "")
                                     if(length(y) == 0){
                                       return(NA)
                                     }
                                     y
                                     
                                   })
                     x <- stringr::str_replace(string = x, 
                                               pattern = c("CAS.ID\\:", "HMDB.ID\\:", "KEGG.ID\\:", "Lab.ID\\:",
                                                           "Adduct\\:", "mz.error\\:", "mz.match.score\\:",
                                                           "RT.error\\:", "RT.match.score\\:",
                                                           "CE\\:", "SS\\:", "Total.score\\:", "Database\\:"), 
                                               replacement = "")
                     x <- matrix(c(Compound.name, x), nrow = 1, byrow = TRUE)
                     colnames(x) <- c("Compound.name", "CAS.ID", "HMDB.ID", "KEGG.ID", "Lab.ID",
                                      "Adduct", "mz.error", "mz.match.score",
                                      "RT.error", "RT.match.score",
                                      "CE", "SS", "Total.score", "Database")
                     return(as.data.frame(x, stringsAsFactors = FALSE))
                   }
                 })
               
               identification <- do.call(rbind, identification)
               identification.table <- tibble::as_tibble(identification.table)
               identification.table <- dplyr::select(.data = identification.table, -(Identification))
               identification.table <- cbind(identification.table, identification)
               
               invisible(identification.table)
               # identification1 <- dplyr::bind_rows(identification)
             }
           })