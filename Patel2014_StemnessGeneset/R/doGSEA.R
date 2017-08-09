
doGSEA <- function (GenesSelected, GenesBackground, ENSG2category, k = 3) 
{
  Genes <- unique(c(GenesSelected, GenesBackground))
  ENSG2category <- ENSG2category[ENSG2category[[2]] %in% Genes, 
                                 ]
  I <- which(ENSG2category[[2]] %in% GenesSelected)
  Cats = sprintf("%s___%s", ENSG2category[, 1], ENSG2category[, 3])
  t = table(Cats)
  useCats = names(t[t >= k])
  ENSG2category = ENSG2category[Cats %in% useCats, ]
  Cats = Cats[Cats %in% useCats]
  ENSG2category$id = factor(Cats)
  ENSG2category = unique(ENSG2category)
  T1 = table(ENSG2category$id)
  T2 = table(ENSG2category$id[ENSG2category$gene_id %in% GenesSelected])
  T <- matrix(0, nr = length(T1), nc = 4)
  T[, 1] <- length(Genes)
  T[, 2] <- T1[]
  T[, 3] <- length(GenesSelected)
  T[, 4] <- T2[]
  row.names(T) <- names(T1)
  T2 = T
  T2[, 1] = T2[, 1] - T2[, 2]
  T2[, 3] = T2[, 3] - T2[, 4]
  T2 = t(T2)
  PV = apply(T2, 2, function(x) {
    unlist(fisher.test(matrix(x, nr = 2, nc = 2), alternative = "greater")[c("p.value", "estimate")])
  })
  Category = unique(ENSG2category[, c("source", "Category", "Name", "id")])
  colnames(T) = c("Tot", "Ann", "nCL", "nENR")
  m = match(row.names(T), Category$id)
  T = data.frame(
    source = Category$source[m], 
    ID = Category$Category[m], 
    T, 
    oddsRatio = PV[2, ], 
    pval = PV[1, ], 
    padj = p.adjust(PV[1, ], method = "BH"), 
    Name = Category$Name[m], 
    stringsAsFactors = FALSE
  )
  T = T[order(T$pval), ]
  Cat2Gene = tapply(ENSG2category$gene_id, ENSG2category$id, 
                    unique)
  Gene2Cat = tapply(as.character(ENSG2category$id), ENSG2category$gene_id, 
                    unique)
  GSEA <- list(T = T, GenesSelected = GenesSelected, GenesBackground = GenesBackground, 
               Cat2Gene = Cat2Gene, Gene2Cat = Gene2Cat)
  GSEA
}
plotGSEA <- function(GSEA, main = NULL, ID = NULL, maxp = 0.05,
                     nrOfSets = 20, xlim,xaxp=NULL,cex=1.2) {
  I = which((GSEA$T$padj <= maxp) & (log2(GSEA$T$odds) > 0.0))
  if (!is.null(ID)) {
    I <- I[as.character(GSEA$T$ID[I]) %in% ID]
  }
  I <- I[order(GSEA$T$oddsRatio[I])]
  J = which((GSEA$T$padj <= maxp) & (log2(GSEA$T$odds) < 0.0))
  if (!is.null(ID)) {
    J <- J[GSEA$T$ID[J] %in% ID]
  }
  J <- J[order(GSEA$T$oddsRatio[J])]
  I <- c(rev(J),I)
  if (nrOfSets > 0) {
    if (length(I) < nrOfSets) { I <- c(rep(NA, nrOfSets-length(I)),I) }
    if (length(I) > nrOfSets) { I <- I[1:nrOfSets] }
  }
  if (length(I) > 0) {
    name = GSEA$T$Name[I]
    odds = log2(GSEA$T$odds[I])
    odds[is.na(I)] = 0.0
    par(mar=c(5,25,4,1)+0.1,xpd=NA)
    xlim = range(odds,finite=TRUE)
    if (xlim[1] > 0.0) { xlim[1] = 0.0 }
    if (xlim[2] < 0.0) { xlim[2] = 0.0 }
    d = diff(xlim)
    xlim = c(xlim[1]-0.05*d,xlim[2]+0.05*d)
    odds[odds == -Inf] = xlim[1]+0.05*d
    odds[odds == Inf] = xlim[2]-0.05*d
    bp = barplot(odds, horiz=TRUE,xaxp=xaxp, xlim=xlim,
                 col = ifelse(odds < 0,brewer.pal(3,"Pastel1")[1],brewer.pal(3,"Pastel1")[2]),
                 main=main, xlab="log2 odds ratio",cex.axis=cex,cex.lab=cex,cex.main=cex)
    text(rep(0.0,length(name)),bp[,1],name,pos=ifelse(odds < 0,4,2),cex=cex)
  }
  invisible(length(I))
}

hwriteGSEA <- function (GSEA, dir = "GSEA", GeneAnnotation = NULL, linktag = colnames(GeneAnnotation), 
                        padj = 0.05) 
{
  GSEA$T = GSEA$T[GSEA$T$padj <= 0.5, ]
  if (nrow(GSEA$T) == 0) {
    page = openPage(file.path(dir, "allcategories.html"), link.css = "hwriter.css", 
                    body.attributes = list(style = "font-size:60%"))
    hwrite("No categories enriched", heading = 3, page = page)
    closePage(page, splash=FALSE)    
    return()
  }
  Categories = row.names(GSEA$T)
  dir.create(dir, showWarnings = FALSE)
  IDX = c("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">", 
          "<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'><head>", 
          "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></meta>", 
          "<title>GSEA</title>", "</head>", "<frameset cols=\"*,*\">", 
          " <frameset rows=\"*,*\">", "   <frame src=\"allcategories.html\" name=\"allcategories\">", 
          "   <frame src=\"allgenes.html\" name=\"allgenes\">", 
          " </frameset>", " <frameset rows=\"*,*\">", "   <frame src=\"category1.html\" name=\"category\">", 
          "   <frame src=\"gene1.html\" name=\"gene\">", " </frameset>", 
          "</frameset>", "</html>")
  writeLines(IDX, sprintf("%s/index.html", dir))
  T = GSEA$T
  cat2link = sprintf("category%d.html", 1:nrow(T))
  names(cat2link) = row.names(T)
  cat2link[!(names(cat2link) %in% Categories)] = "blank.html"
  page = openPage(file.path(dir, "blank.html"), link.css = "hwriter.css", 
                  body.attributes = list(style = "font-size:60%"))
  closePage(page, splash = FALSE)
  I = order(T$pval)
  T = T[I, ]
  cat2link = cat2link[I]
  row.names(T) <- 1:nrow(T)
  row.bgcolor = rep("", nrow(T))
  row.bgcolor[T$padj <= padj] = "#FFFFB7"
  T$oddsRatio <- sprintf("%0.2f", T$oddsRatio)
  T$pval <- sprintf("%0.2e", T$pval)
  T$padj <- sprintf("%0.2e", T$padj)
  T$ID <- hmakeTag("a", T$ID, href = cat2link, target = "category")
  T$Name <- hmakeTag("a", T$Name, href = cat2link, target = "category")
  file.copy(system.file("images/hwriter.css", package = "hwriter"), 
            file.path(dir, "hwriter.css"))
  page = openPage(file.path(dir, "allcategories.html"), link.css = "hwriter.css", 
                  body.attributes = list(style = "font-size:60%"))
  hwrite("All enriched categories", heading = 3, page = page)
  hwrite(T, page = page, nowrap = 1, row.bgcolor = row.bgcolor)
  closePage(page, splash = FALSE)
  gene = sort(unique(c(GSEA$GenesSelected, GSEA$GenesBackground)))
  row.bgcolor = rep("", length(gene))
  row.bgcolor[gene %in% GSEA$GenesSelected] = "#FFFFB7"
  gene2link = sprintf("gene%d.html", 1:length(gene))
  names(gene2link) = gene
  gene2link[!(names(gene2link) %in% GSEA$GenesSelected)] = "blank.html"
  if (is.null(GeneAnnotation)) {
    gene = matrix(gene, nr = length(gene), nc = 1)
    gene[, 1] <- hmakeTag("a", gene[, 1], href = gene2link, 
                          target = "gene")
  }
  else {
    I = which(!(gene %in% row.names(GeneAnnotation)))
    rn = gene[I]
    gene = GeneAnnotation[gene, ]
    row.names(gene)[I] = rn
    I <- order(gene[, 1])
    gene = gene[I, ]
    row.bgcolor = row.bgcolor[I]
    gene2link = gene2link[I]
    for (t in linktag) {
      gene[[t]] <- hmakeTag("a", gene[[t]], href = gene2link, 
                            target = "gene")
    }
  }
  page = openPage(file.path(dir, "allgenes.html"), link.css = "hwriter.css", 
                  body.attributes = list(style = "font-size:60%"))
  hwrite("All genes", heading = 3, page = page)
  hwrite(gene, page = page, nowrap = 1, row.bgcolor = row.bgcolor)
  closePage(page, splash = FALSE)
  for (cat in Categories) {
    page = openPage(file.path(dir, cat2link[cat]), link.css = "hwriter.css", 
                    body.attributes = list(style = "font-size:60%"))
    hwrite(sprintf("Category %s | %s | %s", GSEA$T[cat, c("source")], 
                   GSEA$T[cat, c("ID")], GSEA$T[cat, c("Name")]), heading = 3, 
           page = page)
    hwrite(sprintf("Tot: %d ; Ann: %d ; nCL: %d ; nENR: %d ; oddsRatio: %0.2f ; pval: %0.2e ; padj=%0.2e", 
                   GSEA$T[cat, "Tot"], GSEA$T[cat, "Ann"], GSEA$T[cat, 
                                                                  "nCL"], GSEA$T[cat, "nENR"], GSEA$T[cat, "oddsRatio"], 
                   GSEA$T[cat, "pval"], GSEA$T[cat, "padj"]), page = page)
    Genes2 = GSEA$Cat2Gene[[cat]]
    row.bgcolor = rep("", length(Genes2))
    row.bgcolor[Genes2 %in% GSEA$GenesSelected] = "#FFFFB7"
    if (is.null(GeneAnnotation)) {
      Genes2 = matrix(Genes2, nr = length(gene), nc = 1)
    }
    else {
      I = which(!(Genes2 %in% row.names(GeneAnnotation)))
      rn = Genes2[I]
      Genes2 = GeneAnnotation[Genes2, ]
      row.names(Genes2)[I] = rn
      I <- order(Genes2[, 1])
      row.bgcolor = row.bgcolor[I]
      Genes2 = Genes2[I, ]
      for (t in linktag) {
        Genes2[[t]] <- hmakeTag("a", Genes2[[t]], href = gene2link[row.names(Genes2)], 
                                target = "gene")
      }
    }
    hwrite(Genes2, page = page, row.bgcolor = row.bgcolor)
    closePage(page, splash = FALSE)
  }
  for (g in GSEA$GenesSelected) {
    page = openPage(file.path(dir, gene2link[g]), link.css = "hwriter.css", 
                    body.attributes = list(style = "font-size:60%"))
    hwrite(sprintf("Gene %s", g), heading = 3, page = page)
    if (!is.null(GeneAnnotation)) {
      hwrite(paste(colnames(GeneAnnotation), GeneAnnotation[g, 
                                                            ], sep = ": ", collapse = " ; "), page = page)
    }
    Cats2 = GSEA$Gene2Cat[[g]]
    Cats2 = Cats2[Cats2 %in% row.names(GSEA$T)]
    T = GSEA$T[Cats2, ]
    if (nrow(T) > 0) {
      I = order(T$pval)
      T = T[I, ]
      row.names(T) <- 1:nrow(T)
      row.bgcolor = rep("", nrow(T))
      row.bgcolor[T$padj <= padj] = "#FFFFB7"
      T$oddsRatio <- sprintf("%0.2f", T$oddsRatio)
      T$pval <- sprintf("%0.2e", T$pval)
      T$padj <- sprintf("%0.2e", T$padj)
      T$ID <- hmakeTag("a", T$ID, href = cat2link[row.names(T)], 
                       target = "category")
      T$Name <- hmakeTag("a", T$Name, href = cat2link[row.names(T)], 
                         target = "category")
      row.names(T) <- 1:nrow(T)
    }
    hwrite(T, page = page)
    closePage(page, splash = FALSE)
  }
}

hwriteGSEAtext <- function (GSEA, file, GeneAnnotation = NULL, linktag = colnames(GeneAnnotation), 
                            padj = 0.5) 
{
  T = GSEA$T
  T = T[GSEA$T$padj <= padj, ]
  if (nrow(T) > 0) {
    I = order(T$pval)
    T = T[I, ]
    row.names(T) <- 1:nrow(T)
    T$oddsRatio <- sprintf("%0.2f", T$oddsRatio)
    T$pval <- sprintf("%0.2e", T$pval)
    T$padj <- sprintf("%0.2e", T$padj)
    T$ID <- T$ID
    T$Name <- T$Name
    write.table(T, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

tableGSEA <- function (GSEA, file, GeneAnnotation = NULL, linktag = colnames(GeneAnnotation), 
                            padj = 0.5) 
{
  T = GSEA$T
  T = T[GSEA$T$padj <= padj, ]
  if (nrow(T) > 0) {
    I = order(T$pval)
    T = T[I, ]
    row.names(T) <- 1:nrow(T)
    T$oddsRatio <- sprintf("%0.2f", T$oddsRatio)
    T$pval <- sprintf("%0.2e", T$pval)
    T$padj <- sprintf("%0.2e", T$padj)
    T$ID <- T$ID
    T$Name <- T$Name
    T
  }
}

hwriteGSEArevigo <- function (GSEA, file, GeneAnnotation = NULL, linktag = colnames(GeneAnnotation), 
                              padj = 0.5) 
{
  T = GSEA$T
  T = T[T$padj <= padj, ]
  if (nrow(T) > 0) {
    T = T[T$source %in% c("GO.BP", "GO.CC", "GO.MF"), ]
    T = T[T$oddsRatio > 1, ]
    I = order(T$pval)
    T = T[I, ]
    row.names(T) <- 1:nrow(T)
    T$pval <- sprintf("%0.2e", T$pval)
    T$ID <- T$ID
    T = T[, c("ID", "pval")]
    write.table(T, file = file, sep = "\t", quote = FALSE, row.names = FALSE, 
                col.names = FALSE)
  }
}

hwriteGSEAsimple <- function(GSEA, dir = "GSEA", GeneAnnotation=NULL, linktag=colnames(GeneAnnotation), padj=0.05) {
  GSEA$T = GSEA$T[GSEA$T$padj <= 0.5,]
  Categories = row.names(GSEA$T)
  
  dir.create(dir, showWarnings=FALSE)
  
  ##################
  ## all enriched categories
  T = GSEA$T
  if (nrow(T) > 0) {
    #   cat2link = sprintf("category%d.html", 1:nrow(T))
#   names(cat2link) = row.names(T)
#   cat2link[!(names(cat2link) %in% Categories)] = "blank.html"
#   page = openPage(file.path(dir, "index.html"), link.css="hwriter.css",body.attributes=list(style="font-size:60%"))
#   closePage(page, splash=FALSE)
  I = order(T$pval)
  T = T[I,]
#   cat2link = cat2link[I]
  row.names(T) <- 1:nrow(T)
  row.bgcolor=rep("", nrow(T))
  row.bgcolor[T$padj <= padj] = "#FFFFB7"
  T$oddsRatio <- sprintf("%0.2f", T$oddsRatio)
  T$pval <- sprintf("%0.2e", T$pval)
  T$padj <- sprintf("%0.2e", T$padj)
  T$ID <- T$ID
  T$Name <- T$Name
  file.copy(system.file("images/hwriter.css",package="hwriter"),file.path(dir,"hwriter.css"))
  page = openPage(file.path(dir, "index.html"), link.css="hwriter.css",body.attributes=list(style="font-size:60%"))
  hwrite("All enriched categories", heading=3, page=page)
  ## row.bgcolor=list('#aaffaa', '3'='#ffffaa', '5'=colors)
  hwrite(T, page=page, nowrap=1, row.bgcolor=row.bgcolor)
  closePage(page, splash=FALSE)
  }
}

