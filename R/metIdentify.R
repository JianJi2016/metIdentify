#' @title metIdentify
#' @description Identify peaks based on MS/MS database.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param ms1.data The name of ms1 peak table (csv format). Column 1 is "name", Column 2 is
#' "mz" and column is "rt" (second).
#' @param ms2.data MS2 data, they must be mgf of msp files. For example, ms2.data = c("test.mgf", "test2.msp").
#' @param ms1.ms2.match.mz.tol MS1 peak and MS2 spectrum matching m/z tolerance.
#' @param ms1.ms2.match.rt.tol MS1 peak and MS2 spectrum matching RT tolerance.
#' @param ms1.match.ppm Precursor match ppm tolerance.
#' @param ms2.match.ppm Fragment ion match ppm tolerance.
#' @param mz.ppm.thr m/z
#' @param ms2.match.tol MS2 match tolerance.
#' @param fraction.weight Fraction weight.
#' @param dp.forward.weight Forward DP weight.
#' @param dp.reverse.weight Reverse DP weight.
#' @param rt.match.tol RT match tolerance.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ce Collision energy.
#' @param column "hilic" or "rp".
#' @param ms1.match.weight MS1 match weight.
#' @param rt.match.weight RT match weight.
#' @param ms2.match.weight MS2 match weight.
#' @param path Work directory.
#' @param total.score.tol Total score tolerance.
#' @param candidate.num The number of candidate.
#' @param database MS2 database name.
#' @param threads Number of threads
#' @return A metIdentifyClass object.
#' @export

# setwd("D:/Test/metTools/metIdentify")
# ms1.data <- "peak.table.pos.csv"
# file <- dir()
# ms2.data <- grep("mgf", file, value = TRUE)
# ms1.ms2.match.mz.tol = 25
# ms1.ms2.match.rt.tol = 10
# ms1.match.ppm = 25
# ms2.match.ppm = 30
# mz.ppm.thr = 400
# ms2.match.tol = 0.3
# rt.match.tol = 300000000
# polarity = "positive"
# column = "rp"
# ms1.match.weight = 0.3
# rt.match.weight = 0
# ms2.match.weight = 0.7
# total.score.tol = 0.3
# candidate.num = 3
# path = "."
# database = "msDatabase0.0.1"
# threads = 6
# ce = "all"
# # 
# object.ms <- metIdentify(ms1.data = ms1.data,
#                                ms2.data = ms2.data,
#                            ms2.match.tol = ms2.match.tol,
#                          ms1.match.weight = ms1.match.weight,
#                          rt.match.weight = rt.match.weight,
#                          ms2.match.weight = ms2.match.weight,
#                          rt.match.tol = rt.match.tol, 
#             database = database,
#             ce = ce,
#             total.score.tol = total.score.tol,
#             threads = threads)

setGeneric(name = "metIdentify",
           def = function(ms1.data, ##csv format
                          ms2.data,##only msp and mgf and mz(X)ML are supported
                          ms1.ms2.match.mz.tol = 25,
                          ms1.ms2.match.rt.tol = 10,
                          ms1.match.ppm = 25,
                          ms2.match.ppm = 30,
                          mz.ppm.thr = 400,
                          ms2.match.tol = 0.5,
                          fraction.weight = 0.3,
                          dp.forward.weight = 0.6,
                          dp.reverse.weight = 0.1,
                          rt.match.tol = 30,
                          polarity = c("positive", "negative"),
                          ce = "30",
                          column = c("hilic", "rp"),
                          ms1.match.weight = 0.25,
                          rt.match.weight = 0.25,
                          ms2.match.weight = 0.5,
                          path = ".",
                          total.score.tol = 0.5,
                          candidate.num = 3,
                          database,
                          threads = 3){
             ###Check data
             if(missing(database)){
               stop("No database is provided.\n")
             }
             
             ##parameter specification
             polarity <- match.arg(polarity)
             column <- match.arg(column)
              ##check ms1.file and ms2.file
             file <- dir(path)
             
             if(!all(ms1.data %in% file)) {
               stop("MS1 data is not in the directory, please check it.\n") 
             }
             
             if(!all(ms2.data %in% file)) {
               stop("Some MS2 data are not in the directory, please check it.\n") 
             }
             
             if(!all(database %in% file)) {
               stop("Database is not in this directory, please check it.\n") 
             }
             
             #load MS2 database
             database.name <- database
             load(file.path(path, database.name))
             database <- get(database.name)
             
             ce.list.pos <- unique(unlist(lapply(database@spectra.data$Spectra.positive, names)))
             ce.list.neg <- unique(unlist(lapply(database@spectra.data$Spectra.negative, names)))
             ce.list <- ifelse(polarity == "positive", ce.list.pos, ce.list.neg)
             if(all(ce %in% ce.list) & ce != "all"){
               stop("All ce values you set are not in database. Please check it.\n")
               ce <- ce[ce %in% ce.list]
             }
             rm(list = c("ce.list.pos", "ce.list.neg", "ce.list"))
             
             ##ce values
             if(all(ce != "all")){
               if(polarity == "positive"){
                 ce.list <- unique(unlist(lapply(database@spectra.data$Spectra.positive, function(x){
                   names(x)
                 })))
                 if(length(grep("Unknown", ce.list)) > 0){
                   ce <- unique(c(ce, grep(pattern = "Unknown", ce.list, value = TRUE)))
                 }
               }else{
                 ce.list <- unique(unlist(lapply(database@spectra.data$Spectra.negative, function(x){
                   names(x)
                 })))
                 if(length(grep("Unknown", ce.list)) > 0){
                   ce <- unique(c(ce, grep(pattern = "Unknown", ce.list, value = TRUE)))
                 }
               }               
             }
             
             ##RT in database or not
             if(!database@database.info$RT){
               cat("No RT information in database. The weight of RT have been set as 0.\n")
             }
             #------------------------------------------------------------------
             ##load adduct table
             if(polarity == "positive" & column == "hilic"){
               data("hilic.pos", envir = environment())
               adduct.table <- hilic.pos
             }
             
             if(polarity == "positive" & column == "rp"){
               data("rp.pos", envir = environment())
               adduct.table <- rp.pos
             }
             
             if(polarity == "negative" & column == "hilic"){
               data("hilic.neg", envir = environment())
               adduct.table <- hilic.neg
             }
             
             if(polarity == "negative" & column == "rp"){
               data("rp.neg", envir = environment())
               adduct.table <- rp.neg
             }
             
             if(all(c("ms1.info", "ms2.info") %in% file)){
               cat("Use old data\n")
               load(file.path(path, "ms1.info"))
               load(file.path(path, "ms2.info"))
             }else{
               ##read MS2 data
               cat("Reading MS2 data...\n")
               ms2.data.name <- ms2.data
               temp.ms2.type <- stringr::str_split(string = ms2.data.name, 
                                                   pattern = "\\.")[[1]]
               temp.ms2.type <- temp.ms2.type[length(temp.ms2.type)]
               
               if(temp.ms2.type %in% c("mzXML", "mzML")){
                 ms2.data <- readMZXML(file = ms2.data.name, threads = threads)
               }else{
                 ms2.data <- lapply(ms2.data.name, function(temp.ms2.data){
                   temp.ms2.type <- stringr::str_split(string = temp.ms2.data, 
                                                       pattern = "\\.")[[1]]
                   temp.ms2.type <- temp.ms2.type[length(temp.ms2.type)]
                   if(!temp.ms2.type %in% c("mgf", "msp")) stop("We only support mgf or msp.\n")
                   if(temp.ms2.type == "msp"){
                     temp.ms2.data <- readMSP(file = temp.ms2.data)
                   }else{
                     temp.ms2.data <- readMGF(file = temp.ms2.data)
                   }
                   temp.ms2.data
                 })
                 
                 names(ms2.data) <- ms2.data.name
                 ###prepare data for metIdentification function
                 cat("Preparing MS2 data for identification...\n")
                 ms2.data <- mapply(FUN = function(temp.ms2.data, temp.ms2.data.name){
                   temp.ms2.data <- lapply(temp.ms2.data, function(x){
                     info <- x$info 
                     info <- data.frame(name = paste("mz", info[1],"rt", info[2], sep = ""),
                                        "mz" = info[1], "rt" = info[2], 
                                        "file" = temp.ms2.data.name,
                                        stringsAsFactors = FALSE)
                     rownames(info) <- NULL
                     x$info <- info
                     x
                   })
                   temp.ms2.data
                 },
                 temp.ms2.data = ms2.data,
                 temp.ms2.data.name = ms2.data.name)
                 
                 if(class(ms2.data) == "matrix"){
                   ms2.data <- ms2.data[,1]
                 }else{
                   ms2.data <- do.call(what = c, args = ms2.data)   
                 }  
               }
               
               ms1.info <- lapply(ms2.data, function(x){
                 x[[1]]
               })  
               
               ms2.info <- lapply(ms2.data, function(x){
                 x[[2]]
               })
               
               ms1.info <- do.call(what = rbind, args = ms1.info)
               ms1.info <- as.data.frame(ms1.info)
               rownames(ms1.info) <- NULL
               
               duplicated.name <- unique(ms1.info$name[duplicated(ms1.info$name)])
               if(length(duplicated.name) > 0){
                 lapply(duplicated.name, function(x){
                   ms1.info$name[which(ms1.info$name == x)] <- paste(x, c(1:sum(ms1.info$name == x)), sep = "_")
                 })
               }
               
               names(ms2.info) <- ms1.info$name
               ##save intermediate data
               save(ms1.info, file = file.path(path, "ms1.info"), compress = "xz")
               save(ms2.info, file = file.path(path, "ms2.info"), compress = "xz")
             }
             
             
             if(!missing(ms1.data)){
               cat("Matching peak table with MS2 spectrum...\n")
               ms1.data <- readr::read_csv(file = ms1.data, 
                                           col_types = readr::cols())
               colnames(ms1.data)[1:3] <- c("name", "mz", "rt")
               match.result <- SXTMTmatch(data1 = ms1.data[,c(2,3)], 
                                          data2 = ms1.info[,c(2,3)], 
                                          mz.tol = ms1.ms2.match.mz.tol, 
                                          rt.tol = ms1.ms2.match.rt.tol, 
                                          rt.error.type = "abs")
               if(is.null(match.result)) return("No peaks are matched with MS2 spectra.\n")
               if(nrow(match.result) == 0) return("No peaks are matched with MS2 spectra.\n")
               cat(length(unique(match.result[,1])),"out of", nrow(ms1.data), "peaks have MS2 spectra.\n")
               
               ###if one peak matches multiple peaks, select the more relibale MS2 spectrum
               cat("Selecting the most intense MS2 spectrum for each peak...\n")
               
               temp.idx <- unique(match.result[,1])
               
               match.result <- lapply(temp.idx, function(idx){
                 idx2 <- match.result[which(match.result[,1] == idx), 2]
                 if(length(idx2) == 1){
                   return(c(idx, idx2))
                 }else{
                   temp.ms2.info <- ms2.info[idx2]
                   return(c(idx, idx2[which.max(unlist(lapply(temp.ms2.info, function(y){
                     y <- y[order(y[,2], decreasing = TRUE),,drop = FALSE]
                     if(nrow(y) > 5) y <- y[1:5,]
                     sum(y[,2])
                   })))]
                   )
                   )
                 }
               })
               
               match.result <- do.call(rbind, match.result)
               match.result <- as.data.frame(match.result)
               match.result <- data.frame(match.result, 
                                          ms1.data$name[match.result$V1],
               ms1.info$name[match.result$Index2], stringsAsFactors = FALSE)
               colnames(match.result) <- c("Index1.ms1.data", "Index.ms2.spectra", 
                                           "MS1.peak.name", "MS2.spectra.name")
               ms1.info <- ms1.info[unique(match.result[,2]), , drop = FALSE]
               ms2.info <- ms2.info[unique(match.result[,2])]
               
               match.result$Index.ms2.spectra <- match(match.result$MS2.spectra.name, ms1.info$name)
               save(match.result, file = file.path(path, "match.result"), compress = "xz")
             }else{
               stop("Please provide MS1 data name.\n")
             }
             
             ms2Matchresult <- metIdentification(ms1.info = ms1.info, 
                                                 ms2.info = ms2.info, 
                                                 polarity = polarity, 
                                                 ce = ce, 
                                                 database = database, 
                                                 ms1.match.ppm = ms1.match.ppm, 
                                                 ms2.match.ppm = ms2.match.ppm,
                                                 mz.ppm.thr = mz.ppm.thr,
                                                 ms2.match.tol = ms2.match.tol, 
                                                 rt.match.tol = rt.match.tol, 
                                                 column = column, 
                                                 ms1.match.weight = ms1.match.weight, 
                                                 rt.match.weight = rt.match.weight,
                                                 ms2.match.weight = ms2.match.weight,
                                                 total.score.tol = total.score.tol, 
                                                 candidate.num = candidate.num,
                                                 adduct.table = adduct.table,
                                                 threads = threads,
                                                 fraction.weight = fraction.weight,
                                                 dp.forward.weight = dp.forward.weight,
                                                 dp.reverse.weight = dp.reverse.weight
                                                 )
             
             return.result <- new(Class = "metIdentifyClass",
                                  ms1.data = ms1.data,
                                  ms1.info = ms1.info,
                                  ms2.info = ms2.info,
                                  identification.result = ms2Matchresult,
                                  match.result = match.result,
                                  adduct.table = adduct.table,
                                  ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
                                  ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
                                  ms1.match.ppm = ms1.match.ppm,
                                  ms2.match.ppm = ms2.match.ppm,
                                  ms2.match.tol = ms2.match.tol,
                                  rt.match.tol = rt.match.tol,
                                  polarity = polarity,
                                  ce = paste(ce, collapse = ";"),
                                  column = column,
                                  ms1.match.weight = ms1.match.weight,
                                  rt.match.weight = rt.match.weight,
                                  ms2.match.weight = ms2.match.weight,
                                  path = path,
                                  total.score.tol = total.score.tol,
                                  candidate.num = candidate.num,
                                  database = database.name,
                                  threads = threads)
             cat("All is done.\n")
             return(return.result)
           })

###S4 class for function metIdentification
setClass(Class = "metIdentifyClass", 
         representation(ms1.data = "data.frame",
                        ms1.info = "data.frame",
                        ms2.info = "list",
                        identification.result = "list",
                        match.result = "data.frame",
                        adduct.table = "data.frame",
                        ms1.ms2.match.mz.tol = "numeric",
                        ms1.ms2.match.rt.tol = "numeric",
                        ms1.match.ppm = "numeric",
                        ms2.match.ppm = "numeric",
                        ms2.match.tol = "numeric",
                        rt.match.tol = "numeric",
                        polarity = "character",
                        ce = "character",
                        column = "character",
                        ms1.match.weight = "numeric",
                        rt.match.weight = "numeric",
                        ms2.match.weight = "numeric",
                        path = "character",
                        total.score.tol = "numeric",
                        candidate.num = "numeric",
                        database = "character",
                        threads = "numeric")
)


setMethod(f = "show",
          signature = "metIdentifyClass",
          definition = function(object){
            cat("-----------Identifications------------\n")
            cat("(Use getIdentificationTable to get identification table)\n")
            cat("There are", nrow(object@ms1.data), "peaks\n")
            cat(nrow(object@match.result), "peaks have MS2 spectra\n")
            cat("There are",
                length(unique(unlist(lapply(object@identification.result, function(x){
                  x$Compound.name
                })))),
                "metabolites are identified\n"
                )
            if(!is.null(object@identification.result[[1]])){
            cat("There are", length(object@identification.result), "peaks with identification\n")
            }
            
            cat("-----------Parameters------------\n")
            cat("(Use getParam to get all the parameters of this processing)\n")
            cat("Polarity:", object@polarity, "\n")
            cat("Collision energy:", object@ce, "\n")
            cat("database:", object@database, "\n")
            cat("Total score cutoff:", object@total.score.tol, "\n")
            cat("Column:", object@column, "\n")
            cat("Adduct table:\n")
            cat(paste(object@adduct.table$adduct, collapse = ";"))
            # print(head(tibble::as_tibble(object@adduct.table, 5)))
          }
)

#------------------------------------------------------------------------------
#' @title getParams
#' @description Get parameters from a metIdentifyClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @return A data.frame contains all the parameters of this metIdentifiyClass object.
#' @export
setGeneric(name = "getParams", 
           def = function(object){
             if(class(object) != "metIdentifyClass") stop("Only for metIdentifyClass\n")
             data.frame("Parameter" = c("ms1.ms2.match.mz.tol",
                                        "ms1.ms2.match.rt.tol",
                                        "ms1.match.ppm",
                                        "ms2.match.ppm",
                                        "ms2.match.tol",
                                        "rt.match.tol",
                                        "polarity",
                                        "ce",
                                        "column",
                                        "ms1.match.weight",
                                        "rt.match.weight",
                                        "ms2.match.weight",
                                        "path",
                                        "total.score.tol",
                                        "candidate.num",
                                        "database",
                                        "threads"
                                        ),
                       "Meaning" = c("MS1 features & MS spectra matching mz tolerance (ppm)",
                                        "MS1 features & MS spectra matching RT tolerance (s)",
                                        "MS1 match tolerance (ppm)",
                                        "MS2 fragment match tolerance (ppm)",
                                        "MS2 match tolerance",
                                        "RT match tolerance (s)", 
                                        "Polarity", 
                                        "Collision energy",
                                        "Column",
                                        "MS1 match weight",
                                        "RT match weight",
                                        "MS2 match weight",
                                        "Work directory",
                                        "Total score tolerance",
                                        "Candidate number",
                                        "MS2 database", 
                                        "Thread number"),
                         "Value" = c(object@ms1.ms2.match.mz.tol,
                                    object@ms1.ms2.match.rt.tol,
                                    object@ms1.match.ppm,
                                    object@ms2.match.ppm,
                                    object@ms2.match.tol,
                                    object@rt.match.tol,
                                    object@polarity,
                                    object@ce,
                                    object@column,
                                    object@ms1.match.weight,
                                    object@rt.match.weight,
                                    object@ms2.match.weight,
                                    object@path,
                                    object@total.score.tol,
                                    object@candidate.num,
                                    object@database,
                                    object@threads)
                        )
           })



##------------------------------------------------------------------------------
#' @title getIdentificationTable
#' @description Get identification table from a metIdentifyClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @param type The stype of identification table.
#' @return A identification table (data.frame).
#' @export
setGeneric(name = "getIdentificationTable", 
           def = function(object,
                          type = c("old", "new")
                          ){
             if(class(object) != "metIdentifyClass") stop("Only for metIdentifyClass\n")
             type <- match.arg(type)
             database <- object@database
             identification.result <- object@identification.result
             if(is.null(identification.result[[1]])){
               return(NULL)
             }
             ##add database information
             identification.result <- lapply(identification.result, function(x){
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
                                                   "CE", "SC", "Total.score", "Database"))
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
             identification.table
           })

##------------------------------------------------------------------------------
#' @title getIdenInfo
#' @description Get identification information from a metIdentifyClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @param which.peak Peak name or "all".
#' @param database Database used.
#' @return A identification table (data.frame).
#' @export
setGeneric(name = "getIdenInfo",
           def = function(object,
                          which.peak,
                          database){
             if(missing(object) | missing(which.peak) | missing(database)){
               stop("Please provide the object, which.peak and database.\n")
             }
             if(class(object) != "metIdentifyClass") stop("Only for metIdentifyClass\n")
             if(class(database) != "databaseClass") stop("Only for databaseClass\n")

             identification.result <- object@identification.result

             which.peak <- as.character(which.peak)
             if(!which.peak %in% object@ms1.data$name){
               stop(which.peak, " is not in peak table, please check it.\n")
             }
               
             if(is.na(match(which.peak, object@match.result$MS1.peak.name))){
               cat("The peak has no MS2 spectrum.\n")
               return()
             }

             if(is.na(match(object@match.result$MS2.spectra.name[match(which.peak, object@match.result$MS1.peak.name)],
                   names(identification.result)))){
               cat("The peak hsa no identification result.\n")
               return()
             }

             temp <- match(object@match.result$MS2.spectra.name[match(which.peak, object@match.result$MS1.peak.name)],
                   names(identification.result))
             temp <- identification.result[[temp]]
             temp <- data.frame(temp, database@spectra.info[match(temp$Lab.ID, database@spectra.info$Lab.ID),
                                   setdiff(colnames(database@spectra.info), colnames(temp))
                                   ], stringsAsFactors = FALSE)
             temp <- tibble::as_tibble(temp)
             temp
           })

##------------------------------------------------------------------------------
#' @title ms2plot
#' @description Get MS2 match plots from a metIdentifyClass object.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @param database Used database.
#' @param which.peak Peak name(s) or "all".
#' @param ppm.tol MS2 fragment match ppm.
#' @param mz.ppm.thr m/z
#' @param path Work directory.
#' @param width inch.
#' @param height inch.
#' @param interaction.plot TRUE or FALSE.
#' @param range.mz m/z range.
#' @param range.int Relative intensity range.
#' @param xlab Title of x axis.
#' @param ylab Title of y axis.
#' @param col.lib Colour of database MS2 spectrum.
#' @param col.exp Colour of experimental MS2 spectrum.
#' @param title.size Font size of title.
#' @param lab.size Font size of title of axis.
#' @param axis.text.size Font size of axis text.
#' @param legend.title.size Legend title size.
#' @param legend.text.size Legend text size.
#' @param figure.type "pdf" or "png".
#' @param threads threads
#' @param one.folder Output all figure in one folder?
#' @return A or all ms2 match plot(s).
#' @export 


# object = object
# database = nistDatabase0.0.1
# which.peak = "9.25_379.2487n"
# path = "."
# width = 20
# height = 7
# interaction.plot = FALSE
# range.int = c(-1, 1)
# xlab = "Mass to charge ratio (m/z)"
# ylab = "Relative intensity"
# col.lib = "red"
# col.exp = "black"
# col.filtered = "gray"
# title.size = 15
# lab.size = 12
# axis.text.size =12
# legend.title.size = 12
# legend.text.size = 10
# figure.type = "png"
# threads = 3

setGeneric(name = "ms2plot", 
           def = function(object,
                          database,
                          which.peak = "all",
                          ppm.tol = 30,
                          mz.ppm.thr = 400,
                          path = ".",
                          width = 20, 
                          height = 8,
                          interaction.plot = FALSE,
                          range.mz,
                          range.int = c(-1, 1),
                          xlab = "Mass to charge ratio (m/z)",
                          ylab = "Relative intensity",
                          col.lib = "red",
                          col.exp = "black",
                          title.size = 15,
                          lab.size = 12,
                          axis.text.size =12,
                          legend.title.size = 12,
                          legend.text.size = 10,
                          figure.type = c("png", "pdf"),
                          threads = 3,
                          one.folder = TRUE
           ){
             # 
             if(class(object) != "metIdentifyClass") stop("Only for metIdentifyClass\n")
             if(which.peak == "all"){
               which.peak <- object@ms1.data$name
             }
             identification.result <- object@identification.result
             polarity <- object@polarity
             figure.type <- match.arg(figure.type)
             ##-------------------------------------------------------------------
             ##only for one peak
             if(all(which.peak != "all") & length(which.peak) == 1){
               which.peak <- as.character(which.peak)
               if(!which.peak %in% object@ms1.data$name)
                 stop(which.peak, " is not in peak table, please check it.\n")
               ms2.spectra.name <- object@match.result$MS2.spectra.name[match(which.peak, 
                                                                              object@match.result$MS1.peak.name)]
               if(is.na(ms2.spectra.name)){
                 cat(which.peak, "hsa no MS2 spectrum.\n")
                 return()
               }
               temp.idx <- which(names(identification.result) == ms2.spectra.name)
               if(length(temp.idx) == 0){
                 cat(which.peak, "hsa no identification.\n")
                 return()
               }
               matched.info <- identification.result[[temp.idx]]
                 
                 if(nrow(matched.info) > 1){
                   cat("There are", nrow(matched.info), "identifications.\n")
                   cat(paste(paste(c(1:nrow(matched.info)), as.character(matched.info[,1]), sep = ":"), collapse = "\n"))
                   cat("\n")
                   which.identification <- "test"
                   while(is.na(which.identification) | !which.identification %in% c(1:length(matched.info))){
                     which.identification <- readline(prompt = "Which identification (index: number)?")
                     which.identification <- as.numeric(which.identification)
                   }
                   matched.info <- unlist(matched.info[which.identification, , drop = TRUE])
                 }else{
                   matched.info <- unlist(matched.info[1, , drop = TRUE])
                 }
               
                 lib.spectrum <- getMS2spectrum(lab.id = matched.info["Lab.ID"], 
                                                database = database, 
                                                polarity = polarity, 
                                                ce = matched.info["CE"])
                 exp.spectrum <- object@ms2.info[[match(ms2.spectra.name, names(object@ms2.info))]]
                 if(missing(range.mz)){
                   range.mz <- range(c(lib.spectrum[, "mz"], exp.spectrum[, "mz"]))  
                 }
                 
                 plot <- plotMS2match(matched.info = matched.info,
                                           ppm.tol = ppm.tol,
                                           mz.ppm.thr = mz.ppm.thr,
                                           exp.spectrum = exp.spectrum,
                                           lib.spectrum = lib.spectrum,
                                           polarity = polarity,
                                           xlab = xlab,
                                           ylab = ylab,
                                           col.lib = col.lib,
                                           col.exp = col.exp,
                                           ce = matched.info["CE"],
                                           title.size = title.size,
                                           lab.size = lab.size,
                                           axis.text.size = axis.text.size,
                                           legend.title.size = legend.title.size,
                                           legend.text.size = legend.text.size,
                                           database = database
                 )
                 if(interaction.plot){
                   plot <- plotly::ggplotly(p = plot)
                 }
                 plot
             }else{
               ##output all MS2 match
               dir.create(path)
               path <- file.path(path, "ms2_match_plot")
               dir.create(path)
               
               if(all(which.peak != "all")){
                 if(!all(which.peak %in% object@ms1.data$name)){
                   stop("Some peaks are not in MMS1 peak table, please check them.\n")
                 }
                 ms2.spectra.name <- object@match.result$MS2.spectra.name[match(which.peak, object@match.result$MS1.peak.name)]
                 which.peak <- which.peak[!is.na(ms2.spectra.name)]
                 ms2.spectra.name <- ms2.spectra.name[!is.na(ms2.spectra.name)]
                 if(length(ms2.spectra.name) == 0){
                   stop("All peaks have no MS2 spectra.\n")
                 }
                 anno.idx <- match(ms2.spectra.name, names(object@identification.result))  
                 if(all(is.na(anno.idx))){
                   stop("All peaks have no identifications.\n")
                 }
                 
                 which.peak <- which.peak[!is.na(anno.idx)]
                 ms2.spectra.name <- ms2.spectra.name[!is.na(anno.idx)]    
                 anno.idx <- anno.idx[!is.na(anno.idx)]
               }
               cat("There are", length(anno.idx), "peaks with identifications.\n")
               
               if(length(anno.idx) == 0) {
                 return(NULL)
               }
               temp.fun <- function(anno.idx, 
                                    identification.result,
                                    ms2.info,
                                    match.result,
                                    database,
                                    ppm.tol = 30,
                                    mz.ppm.thr = 400,
                                    col.lib = "red",
                                    col.exp = "black",
                                    polarity = c("positive", "nagative"),
                                    range.int = c(-1, 1), 
                                    xlab = "Mass to charge ratio (m/z)",
                                    ylab = "Relative intensity", 
                                    title.size = 15, 
                                    lab.size = 12, 
                                    axis.text.size = 12, 
                                    legend.title.size = 12,  
                                    legend.text.size = 10,
                                    plotMS2match,
                                    getMS2spectrum
               ){
                 matched.info <- identification.result[[anno.idx]]  
                 temp.ms2.spectrum.name <- names(identification.result)[anno.idx]
                 temp.peak.name <- match.result$MS1.peak.name[match(temp.ms2.spectrum.name,
                                                                    match.result$MS2.spectra.name)]
                 
                 if(!one.folder){
                   temp.path <- file.path(path, 
                                          stringr::str_replace_all(string = temp.peak.name, 
                                                                   pattern = "/", replacement = "_"))
                   dir.create(temp.path)
                 }
                
                 matched.info <- apply(matched.info, 1, list)
                 matched.info <- lapply(matched.info, unlist)
                    
                     non.meaning <- lapply(matched.info, function(temp.matched.info){
                       if(one.folder) {
                         temp.file.name <- 
                           file.path(path,
                                    stringr::str_c(
                                      paste(stringr::str_replace_all(string = temp.peak.name, pattern = "/", replacement = "_"),
                                            paste(as.character(temp.matched.info[c("Total.score", "Lab.ID", "Adduct")]), collapse = ";"),
                                            sep = ";"), ".", figure.type ,sep = "")
                         )
                       }else{
                         temp.file.name <- file.path(temp.path,
                                                     stringr::str_c(paste(as.character(temp.matched.info[c("Total.score", "Lab.ID", "Adduct")]), collapse = ";"), 
                                                                    ".", figure.type ,sep = "")
                         )
                       }
                       
                       
                       lib.spectrum <- getMS2spectrum(lab.id = temp.matched.info["Lab.ID"], 
                                                      database = database, 
                                                      polarity = polarity, 
                                                      ce = temp.matched.info["CE"])
                       exp.spectrum <- ms2.info[[match(match.result$MS2.spectra.name[match(temp.peak.name,
                                                                                           match.result$MS1.peak.name)], 
                                                       names(ms2.info))]]
                       range.mz <- range(c(lib.spectrum[, "mz"], exp.spectrum[, "mz"])) 
                       
                       temp.plot <- plotMS2match(matched.info = temp.matched.info,
                                ppm.tol = ppm.tol,
                                mz.ppm.thr = mz.ppm.thr,
                                exp.spectrum = exp.spectrum,
                                lib.spectrum = lib.spectrum,
                                polarity = polarity,
                                xlab = xlab,
                                ylab = ylab,
                                col.lib = col.lib,
                                col.exp = col.exp,
                                ce = temp.matched.info["CE"],
                                title.size = title.size,
                                lab.size = lab.size,
                                axis.text.size = axis.text.size,
                                legend.title.size = legend.title.size,
                                legend.text.size = legend.text.size,
                                database = database
                       )
                       ggplot2::ggsave(filename = temp.file.name,
                                       plot = temp.plot, 
                                       width = width, height = height)
                       
                     })
                   }
                 
               BiocParallel::bplapply(X = anno.idx, 
                                      FUN = temp.fun, 
                                      BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                        progressbar = TRUE),
                                      identification.result = identification.result,
                                      ms2.info = object@ms2.info,
                                      match.result = object@match.result,
                                      database = database,
                                      ppm.tol = ppm.tol,
                                      mz.ppm.thr = mz.ppm.thr,
                                      col.lib = col.lib,
                                      col.exp = col.exp,
                                      polarity = polarity,
                                      xlab = xlab,
                                      ylab = ylab, 
                                      title.size = title.size, 
                                      lab.size = lab.size, 
                                      axis.text.size = axis.text.size, 
                                      legend.title.size = legend.title.size,  
                                      legend.text.size = legend.text.size,
                                      plotMS2match = plotMS2match,
                                      getMS2spectrum = getMS2spectrum
               )
               cat("All is done.\n")
             }
           })

#------------------------------------------------------------------------------
#' @title whichHasIden
#' @description Get the spectra names with identifications.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @return Spectra names with identifications.
#' @export
setGeneric(name = "whichHasIden", 
           def = function(object){
             if(class(object) != "metIdentifyClass") stop("Only for metIdentifyClass\n")
             temp <- object@match.result[match(names(object@identification.result), object@match.result$MS2.spectra.name),c(3,4)]
             rownames(temp) <- NULL
             temp
           })

#------------------------------------------------------------------------------
#' @title filterIden
#' @description Filter identifications according to mz error, RT error, MS similarity and total score.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metIdentifyClass object.
#' @param ms1.match.ppm MS1 match ppm.
#' @param rt.match.tol RT match tolerance.
#' @param ms2.match.tol MS2 match tolerance.
#' @param total.score.tol Total score tolerance.
#' @return An new metIdentifyClass.
#' @export
setGeneric(name = "filterIden",
           def = function(object,
                          ms1.match.ppm = 25,
                          rt.match.tol = 30,
                          ms2.match.tol = 0.5,
                          total.score.tol = 0.5
                          ){
             if(class(object) != "metIdentifyClass") {
               stop("Only for metIdentifyClass\n")
             }
             
             object@ms1.match.ppm <- ms1.match.ppm
             object@rt.match.tol <- rt.match.tol
             object@ms2.match.tol <- ms2.match.tol
             object@total.score.tol <- total.score.tol
             
             identification.result <- object@identification.result
             
             
             identification.result <- lapply(identification.result, function(x){
               RT.error <- x$RT.error
               RT.error[is.na(RT.error)] <- rt.match.tol - 1
               x <- x[which(x$mz.error < ms1.match.ppm & RT.error < rt.match.tol & 
                              x$SC > ms2.match.tol & x$Total.score > total.score.tol),,drop = FALSE]
             })
             
             temp.idx <- which(unlist(lapply(identification.result, function(x){
               nrow(x) != 0
             })))
             
             identification.result <- identification.result[temp.idx]
             object@identification.result <- identification.result
             object
           })

#' #------------------------------------------------------------------------------
#' #' @title combineObject
#' #' @description Combine identification results from different libries.
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@163.com}
#' #' @param ... More than one metIdentifyClass objects.
#' #' @return An identification table.
#' #' @export
#' setGeneric(name = "combineObject", 
#'            def = function(...){
#'              objects <- list(...)
#'              if(any(unique(unlist(lapply(objects, class))) != "metIdentifyClass")){
#'                stop("Only for metIdentifyClass\n")
#'              }
#'              
#'              if(length(objects) == 1){
#'                cat("Only one object, no combiination.\n")
#'                return(objects[[1]])
#'              }
#'              
#'              identification.tables <- lapply(objects, function(x) x@identification.result)
#'              
#'              if(length(unique(unlist(lapply(identification.tables, function(x){
#'                nrow(x)
#'              })))) != 1){
#'                stop("You objects may be from different data sets, please check them.\n")
#'              }
#'              
#'              
#'              temp.hits.forwards <- lapply(identification.tables, function(x){
#'                x$hits.forward
#'              })
#'              
#'              temp.hits.reverses <- lapply(identification.tables, function(x){
#'                x$hits.reverse
#'              })
#'              
#'              temp.database <- sapply(objects, function(x) x@database)
#'              
#'              temp.hits.forwards <- mapply(function(x, y){
#'                y <- paste("database{", y, "}", sep = "")
#'                x <- unlist(lapply(x, function(z){
#'                  if(z == "") return(z)
#'                  z <- stringr::str_split(string = z, pattern = ";")[[1]]
#'                  z <- stringr::str_c(stringr::str_c(z, y, sep = ""), collapse = ";")
#'                  z
#'                }))
#'                list(x)
#'              },
#'              x = temp.hits.forwards,
#'              y = temp.database)
#'              
#'              
#'              temp.hits.forwards <- do.call(what = cbind, args = temp.hits.forwards)
#'              temp.hits.forwards <- apply(temp.hits.forwards, 1, function(x){
#'                if(all(x == "")) return("")
#'                x <- x[x!=""]
#'                x <- stringr::str_c(x, collapse = ";")
#'                x.score <- stringr::str_extract_all(string = x, 
#'                                                    pattern = "score\\{[0-1]{1}\\.?[0-9]*")[[1]]
#'                x.score <- as.numeric(stringr::str_replace_all(string = x.score,
#'                                                               pattern = "score\\{", replacement = ""))
#'                x <- stringr::str_split(string = x, pattern = ";")[[1]]
#'                x <- x[order(x.score, decreasing = TRUE)]
#'                x <- stringr::str_c(x, collapse = ";")
#'                x
#'              })
#'              
#'              temp.hits.forwards <- unlist(temp.hits.forwards)
#'              
#'              
#'              temp.hits.reverses <- mapply(function(x, y){
#'                y <- paste("database{", y, "}", sep = "")
#'                x <- unlist(lapply(x, function(z){
#'                  if(z == "") return(z)
#'                  z <- stringr::str_split(string = z, pattern = ";")[[1]]
#'                  z <- stringr::str_c(stringr::str_c(z, y, sep = ""), collapse = ";")
#'                  z
#'                }))
#'                list(x)
#'              },
#'              x = temp.hits.reverses,
#'              y = temp.database)
#'              
#'              
#'              temp.hits.reverses <- do.call(what = cbind, args = temp.hits.reverses)
#'              temp.hits.reverses <- apply(temp.hits.reverses, 1, function(x){
#'                if(all(x == "")) return("")
#'                x <- x[x!=""]
#'                x <- stringr::str_c(x, collapse = ";")
#'                x.score <- stringr::str_extract_all(string = x, 
#'                                                    pattern = "score\\{[0-1]{1}\\.?[0-9]*")[[1]]
#'                x.score <- as.numeric(stringr::str_replace_all(string = x.score,
#'                                                               pattern = "score\\{", replacement = ""))
#'                x <- stringr::str_split(string = x, pattern = ";")[[1]]
#'                x <- x[order(x.score, decreasing = TRUE)]
#'                x <- stringr::str_c(x, collapse = ";")
#'                x
#'              })
#'              
#'              temp.hits.reverses <- unlist(temp.hits.reverses)
#'              
#'              temp.nhits.forward <- unlist(lapply(temp.hits.forwards, function(x){
#'                if(x == "") return(0)
#'                length(stringr::str_split(x, pattern = ";")[[1]])
#'              }))
#'              
#'              temp.nhits.reverse <- unlist(lapply(temp.hits.reverses, function(x){
#'                if(x == "") return(0)
#'                length(stringr::str_split(x, pattern = ";")[[1]])
#'              }))
#'              
#'              identification.table <- data.frame("name"= identification.tables[[1]]$name,
#'                                             "mz" = identification.tables[[1]]$mzmed,
#'                                             "rt" = identification.tables[[1]]$rt,
#'                                             "nhits.forward" = temp.nhits.forward,
#'                                             "hits.forward" = temp.hits.forwards,
#'                                             "nhits.reverse" = temp.nhits.reverse,
#'                                             "hits.reverse" = temp.hits.reverses,
#'                                             stringsAsFactors = FALSE
#'              )
#'              identification.table
#'              
#'            })

#' #------------------------------------------------------------------------------
#' #' @title outputIdenTable
#' #' @description Output identification table for each peak.
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@163.com}
#' #' @param ms1.data MS1 peak table, first 3 column must be name, mz and rt.
#' #' @param identification.table Identification table from metIdentifyClass object.
#' #' @param peak.alignment.mz.tol Peaks in ms1 data and MS2 spectra match mz tolerance.
#' #' @return An identification table for each peak.
#' #' @export
#' setGeneric(name = "outputIdenTable", 
#'            def = function(ms1.data,
#'                           identification.table,
#'                           peak.alignment.mz.tol = 25){
#'              ms1.data <- as.data.frame(ms1.data)
#'              identification.table <- as.data.frame(identification.table)
#'              identification.table <- identification.table[which(identification.table$nhits.forward > 0 |
#'                                                           identification.table$nhits.reverse > 0 ), , drop = FALSE]
#'              
#'              if(nrow(identification.table) == 0) {
#'                stop("No identifications in identification.table. Please check it.\n")
#'              }
#'              
#'              match.result <- SXTMTmatch(data1 = ms1.data[,c(2,3,1)], 
#'                                         data2 = identification.table[,c(2,3)], 
#'                                         mz.tol = peak.alignment.mz.tol, 
#'                                         rt.error.type = "abs")
#'              
#'              if(is.null(match.result)){
#'                stop("No peaks in ms1.data are matched with peaks in identification.table.\n")
#'              }
#'              
#'              peak.table.with.identification <- lapply(c(1:nrow(ms1.data)), function(idx){
#'                temp.idx <- which(match.result[,1] == idx)
#'                if(length(temp.idx) == 0) {
#'                  temp <- as.data.frame(matrix(rep(NA, 9), nrow = 1))
#'                  colnames(temp) <- c("matcthed.spectra", "from.which.spectrum", "ms2.spectrm.name",
#'                                      "precursor.mz", "precursor.rt", "nhits.forward", "hits.forward", 
#'                                      "nhits.reverse",
#'                                      "hits.reverse" )
#'                  return(temp)
#'                }
#'                
#'                temp.identification <- identification.table[match.result[temp.idx,2],, drop = FALSE]
#'                temp.identification.spectra.name <- temp.identification[,1]
#'                temp.identification.hits.forward <- temp.identification$hits.forward
#'                temp.identification.hits.reverse <- temp.identification$hits.reverse
#'                if(any(temp.identification.hits.forward != "")){
#'                  temp.identification.forward.score <- stringr::str_extract(string = temp.identification.hits.forward, 
#'                                                                        pattern = "score\\{[0-1]{1}\\.?[0-9]*")
#'                  temp.identification.forward.score <- as.numeric(stringr::str_replace(string = temp.identification.forward.score,
#'                                                                                   pattern = "score\\{", replacement = ""))        
#'                  anno.idx <- which.max(temp.identification.forward.score)
#'                }else{
#'                  temp.identification.reverse.score <- stringr::str_extract(string = temp.identification.hits.reverse, 
#'                                                                        pattern = "score\\{[0-1]{1}\\.?[0-9]*")
#'                  temp.identification.reverse.score <- as.numeric(stringr::str_replace(string = temp.identification.reverse.score,
#'                                                                                   pattern = "score\\{", replacement = ""))        
#'                  anno.idx <- which.max(temp.identification.reverse.score)    
#'                }
#'                
#'                from.which.spectrum <- temp.identification.spectra.name[anno.idx]
#'                
#'                return.result <- temp.identification[anno.idx, , drop = FALSE]
#'                return.result <- data.frame("matcthed.spectra" = stringr::str_c(temp.identification.spectra.name, collapse = ";"),
#'                                            from.which.spectrum,
#'                                            return.result, 
#'                                            stringsAsFactors = FALSE)
#'                colnames(return.result)[c(3:5)] <- c("ms2.spectrm.name", "precursor.mz", "precursor.rt")
#'                return.result
#'              })
#'              
#'              peak.table.with.identification <- do.call(rbind, peak.table.with.identification)
#'              
#'              peak.table.with.identification <- data.frame(ms1.data,peak.table.with.identification, stringsAsFactors = FALSE)
#'              peak.table.with.identification
#'              # write.csv(peak.table.with.identification, 
#'              #           file.path(path, "peak.table.with.identification.csv"), row.names = FALSE)
#'            })



#------------------------------------------------------------------------------
#' @title getMS2spectrum2Object
#' @description Get spectra of peaks.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param peak.name Peak name.
#' @param object metIdentifyClass
#' @return A MS2 spectrum.
#' @export

setGeneric(name = "getMS2spectrum2Object", 
           def = function(peak.name,
                          object){
             if(class(object) != "metIdentifyClass") stop("Only for metIdentifyClass\n")
             if(missing(peak.name)) stop('Please provide peak name.\n')
             object@ms2.info[[which(object@match.result$MS2.spectra.name[match(peak.name, object@match.result$MS1.peak.name)] == names(object@ms2.info))]]
           })













.onAttach <- function(libname, pkgname){
  packageStartupMessage("metIdentify,
More information can be found at https://jaspershen.github.io/metIdentify/
Authors: Xiaotao Shen (shenxt1990@163.com)
Maintainer: Xiaotao Shen.
Version 0.1.1 (20190606)
--------------
o Creatation.")
}

packageStartupMessage("metIdentify,
More information can be found at https://jaspershen.github.io/metIdentify/
Authors: Xiaotao Shen (shenxt1990@163.com)
Maintainer: Xiaotao Shen.
Version 0.1.1 (20190606)
--------------
o Creatation.")