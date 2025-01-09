library(lattice)

manhattan.plot <- function(chr, pos, pvalue, 
                           sig.level = NA, annotate = NULL, ann.default = list(),
                           should.thin = TRUE, thin.pos.places = 2, thin.logp.places = 2, 
                           xlab = "Chromosome", ylab = expression(-log[10](p-value)),
                           col = c("gray", "darkgray"), panel.extra = NULL, pch = 20, cex = 0.8, ...) {
  
  if (length(chr) == 0) stop("chromosome vector is empty")
  if (length(pos) == 0) stop("position vector is empty")
  if (length(pvalue) == 0) stop("pvalue vector is empty")
  
  # Ensure chr is an ordered factor
  if (!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[, drop = TRUE]
  }
  
  # Make sure positions are in kbp
  if (any(pos > 1e6)) pos <- pos / 1e6
  
  # Calculate absolute genomic position from relative chromosomal positions
  posmin <- tapply(pos, chr, min)
  posmax <- tapply(pos, chr, max)
  posshift <- head(c(0, cumsum(posmax)), -1)
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr]
  
  getGenPos <- function(cchr, cpos) {
    p <- posshift[as.character(cchr)] + cpos
    return(p)
  }
  
  # Parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default <- list(x = "peak", y = "peak", adj = NULL, pos = 3, offset = 0.5, 
                        col = NULL, fontface = NULL, fontsize = NULL, show = FALSE) # Label show = FALSE
  parse.label <- function(rawval, groupname) {
    r <- list(text = groupname)
    if (is.logical(rawval)) {
      if (!rawval) { r$show <- FALSE }
    } else if (is.character(rawval) || is.expression(rawval)) {
      if (nchar(rawval) >= 1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if (!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times = length(pvalue)))
  }
  
  ann.settings <- vector("list", length(levels(grp)))
  ann.settings[[1]] <- list(pch = pch, col = col, cex = cex, fill = col, label = label.default)
  
  if (length(ann.settings) > 1) { 
    lcols <- trellis.par.get("superpose.symbol")$col 
    lfills <- trellis.par.get("superpose.symbol")$fill
    for (i in 2:length(levels(grp))) {
      ann.settings[[i]] <- list(pch = pch, 
                                col = lcols[(i-2) %% length(lcols) + 1], 
                                fill = lfills[(i-2) %% length(lfills) + 1], 
                                cex = cex, label = label.default)
      ann.settings[[i]]$label$show <- FALSE # Prevent label display
    }
    names(ann.settings) <- levels(grp)
  }
  
  for (i in 1:length(ann.settings)) {
    if (i > 1) {
      ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)
    }
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  
  if (is.list(annotate) && length(annotate) > 1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if (!is.null(names(annotate[-1])) && all(names(annotate[-1]) != "")) {
      ann.cols <- match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols <- user.cols - 1
    }
    for (i in seq_along(user.cols)) {
      if (!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label <- parse.label(annotate[[user.cols[i]]]$label, 
                                                      levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]] <- modifyList(ann.settings[[ann.cols[i]]], 
                                                annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  # Reduce number of points plotted if should.thin is TRUE
  if (should.thin) {
    thinned <- unique(data.frame(
      logp = round(-log10(pvalue), thin.logp.places), 
      pos = round(genpos, thin.pos.places), 
      chr = chr,
      grp = grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  rm(pos, pvalue)
  gc()
  
  # Custom axis to print chromosome names
  axis.chr <- function(side, ...) {
    if (side == "bottom") {
      panel.axis(side = side, outside = TRUE,
                 at = ((posmax + posmin) / 2 + posshift),
                 labels = levels(chr), 
                 ticks = FALSE, rot = 90,
                 check.overlap = FALSE)
    } else if (side == "top" || side == "right") {
      panel.axis(side = side, draw.labels = FALSE, ticks = FALSE)
    } else {
      axis.default(side = side, ...)
    }
  }
  
  # Make sure the y-lim covers the range
  prepanel.chr <- function(x, y, ...) { 
    A <- list()
    maxy <- ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0))) + .5
    A$ylim <- c(0, maxy)
    A
  }
  
  xyplot(logp ~ genpos, chr = chr, groups = grp,
         axis = axis.chr, ann.settings = ann.settings, 
         prepanel = prepanel.chr, scales = list(axs = "i"),
         panel = function(x, y, ..., getgenpos) {
           if (!is.na(sig.level)) {
             panel.abline(h = -log10(sig.level), lty = 2)
           }
           panel.superpose(x, y, ..., getgenpos = getgenpos)
           if (!is.null(panel.extra)) {
             panel.extra(x, y, getgenpos, ...)
           }
         },
         panel.groups = function(x, y, ..., subscripts, group.number) {
           A <- list(...)
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts]) - 1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts]) - 1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts]) - 1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts]) - 1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
         },
         xlab = xlab, ylab = ylab, 
         panel.extra = panel.extra, getgenpos = getGenPos, ...
  )
}

# Read data from the CSV file
df <- read.csv("./diagenode_v2_results/avpc_vs_group4/data/avpc_vs_group4_DMPs.csv")

# Adjust the factor levels to match the 'chr' prefix format and include 'chrY'
df$chr <- factor(df$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

# Identify significant sites based on the new criteria
up_sites <- which(df$logFC > 1 & df$P.Value < 0.001)
down_sites <- which(df$logFC < -1 & df$P.Value < 0.001)

# Combine both sets of significant sites
significant_sites <- c(up_sites, down_sites)

# Create an annotation vector with default color (gray) for all sites
ann <- rep(1, nrow(df))  # 1 indicates non-significant (default gray)
ann[significant_sites] <- 2  # 2 indicates significant (to be colored differently)

# Convert the annotation vector into a factor
ann <- factor(ann)

# Call the modified manhattan.plot function
manhattan.plot(df$chr, df$pos, df$P.Value, annotate = ann)
