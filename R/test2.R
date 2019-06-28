  # setwd("E:/project/si_wu/1stLayer_annotation_PE/RPLC_pos/in-house/NCE50")
  # ms1.peak.table <- readr::read_csv("PeakTable_pRPLC.csv")
  # ms1.peak.table$`Retention time (min)` <- ms1.peak.table$`Retention time (min)` * 60
  # colnames(ms1.peak.table) <- c("name", "mz", "rt")
  # write.csv(ms1.peak.table, "ms1.peak.table.csv", row.names = FALSE)
  # ms1.data <- "ms1.peak.table.csv"
  # file <- dir()
  # ms2.data <- grep("mgf", file, value = TRUE)
  # ms1.ms2.match.mz.tol = 25
  # ms1.ms2.match.rt.tol = 10
  # ms1.match.ppm = 25
  # ms2.match.ppm = 30
  # mz.ppm.thr = 400
  # ms2.match.tol = 0.4
  # rt.match.tol = 60
  # polarity = "positive"
  # column = "rp"
  # # ms1.match.weight = 0.3
  # # rt.match.weight = 0
  # # ms2.match.weight = 0.7
  # total.score.tol = 0.3
  # candidate.num = 3
  # path = "."
  # database = "msDatabase_rplc0.0.1"
  # threads = 6
  # ce = "all"
  # #
  # result.ms.pos50 <- metIdentify(ms1.data = ms1.data,
  #                                  ms2.data = ms2.data,
  #                                  ms2.match.tol = ms2.match.tol,
  #                                  rt.match.tol = rt.match.tol,
  #                                  database = database,
  #                                  ce = ce,
  #                                  total.score.tol = total.score.tol,
  #                                  threads = threads,
  #                                  column = column)
  # 
  # save(result.ms.pos50, file = "result.ms.pos50")
  # 
  # load(database)
  # ms2plot(object = result.ms.pos50,
  #         database = msDatabase_rplc0.0.1,
  #         which.peak = "all", one.folder = TRUE)
  # 
  # 
  # annotation.table2 <- getIdentificationTable(result.metlin.pos25,
  #                                             result.metlin.pos50,
  #                                             result.nist.pos25,
  #                                             result.nist.pos50,
  #                                             result.hmdb.pos25,
  #                                             result.hmdb.pos50,
  #                                             result.mona.pos25,
  #                                             result.mona.pos50,
  #                                             candidate.num = 3)
  # dim(annotation.table2)
  # annotation.table2 <- annotation.table2[which(!is.na(annotation.table2$Identification)), , drop = FALSE]
  # 
  # annotation.table2 <- annotation.table2[which(!annotation.table2$name %in% annotation.table1$name),]
  # 
  # 
  # write.csv(annotation.table2, file = "annotation.table2.csv", row.names = FALSE)
  # 
  # 
  # 
  # setwd("E:/project/si_wu/1stLayer_annotation_PE/RPLC_pos/ms1")
  # result.mz.pos <- mzIdentify(ms1.data = "ms1.peak.table.csv",
  #                             ms1.match.ppm = 25,
  #                             polarity = "positive", 
  #                             column = "rp", path = ".", 
  #                             candidate.num = 3,
  #                             database = "HMDB.metabolite.data", 
  #                             threads = 5)
  # 
  # 
  # annotation.table3 <- getIdentificationTable2(object = result.mz.pos, 
  #                                              candidate.num = 3)
  # annotation.table3 <- annotation.table3[which(!is.na(annotation.table3$Identification)), , drop = FALSE]
  # annotation.table3 <- annotation.table3[which(!annotation.table3$name %in% unique(c(annotation.table1$name,
  #                                                                             annotation.table2$name))),]
  # 
  # write.csv(annotation.table3, file = "annotation.table3.csv", row.names = FALSE)
  # 
  # 
  # dim(annotation.table1)
  # dim(annotation.table2)
  # dim(annotation.table3)
  # 
  # colnames(annotation.table1)
  # colnames(annotation.table2)
  # colnames(annotation.table3)
  # 
  # annotation.table1 <- data.frame(annotation.table1, Level = 1, stringsAsFactors = FALSE)
  # annotation.table2 <- data.frame(annotation.table2, Level = 2, stringsAsFactors = FALSE)
  # annotation.table3 <- data.frame(annotation.table3, Level = 3, "MS2.spectrum.name" = NA,
  #                                 stringsAsFactors = FALSE)
  # annotation.table3 <- annotation.table3[,c("name", "mz", "rt", "MS2.spectrum.name",
  #                                           "Candidate.number", "Identification", "Level")]
  # 
  # annotation.table <- rbind(annotation.table1, annotation.table2, annotation.table3)
  # 
  # write.csv(annotation.table, file = "annotation.table.csv", row.names = FALSE)






