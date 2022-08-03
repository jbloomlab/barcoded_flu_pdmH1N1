#! /usr/bin/Rscript
### Argument parsing ###########################################################
args <- commandArgs(trailingOnly=TRUE)
csv.filename <- args[1]
plot.filename <- args[2]

### Packages ###################################################################
library(ggplot2)
library(gggenes)
library(ggnewscale)
library(patchwork)

### Helper Fxns ################################################################

# Extract the mutations for all segments for each cell
#
# @param data data.frame containing mutation information
# @param segments vector of segment names to parse
# @param segment.df data.frame containing the coordinates of each segment
# @param pattern prefix pattern for each data.frame column for each segment
#
# @return Returns A data.frame containing mutations, insertions, and deletion
#         locations for each segment.

ParseMutations <- function(data, segments, segment.df, pattern = "mutations_") {
  mutations.df <- list()
  for (i in 1:nrow(x = data)) {
    mutations <- lapply(X = segments, FUN = function(x) {
      ExtractSegmentMutations(
        x = data[i, paste0(pattern, x)], 
        segment = x, 
        segment.df = segment.df, 
        cell = i)
    })
    mutations.df[[i]] <- do.call(what = rbind, args = mutations)
  }
  mutations.df <- do.call(what = rbind, args = mutations.df)
  mutations.df$type <- factor(x = mutations.df$type)
  return(mutations.df)
}

# Extract the mutations on a given segment for a given cell
#
# @param x A vector of mutations
# @param segment Name of the segment mutation is present on
# @param segment.df Data.frame containing the coordinates of each segment
# @param cell Index of cell - added into returned data.frame
#
# @return Returns a data.frame of mutations, insertions, and deletion locations

ExtractSegmentMutations <- function(x, segment, segment.df, cell) {
  if (x == "WT" | x == "Not Detected") {
    return(NULL)
  }
  mutations <- strsplit(x = x, split = " ")[[1]]
  locations <- NULL
  mutation.df <- list()
  deletion.df <- list()
  insertion.df <- list()
  for (mutation in mutations) {
    if (grepl(pattern = "^del", x = mutation)) {
      del.loc <- as.numeric(strsplit(x = gsub(pattern = "del", x = mutation, replacement = ""), split = "to")[[1]])
      del.loc <- del.loc + segment.df[segment.df$segment == segment, "start"]
      deletion.df[[length(deletion.df)+1]] <- data.frame(start = del.loc[1], end = del.loc[2], type = "deletion")
    } else if (grepl(pattern = "^ins", x = mutation)) {
      ins.loc <- gsub(pattern = "ins", x = mutation, replacement = "")
      ins.length <- nchar(gsub(pattern = "([1-9]+)", x = ins.loc, replacement = ""))
      ins.loc <- as.numeric(gsub(pattern = "([A-Z,a-z]+)", x = ins.loc, replacement = "")) + segment.df[segment.df$segment == segment, "start"]
      insertion.df[[length(insertion.df) + 1]] <- data.frame(start = ins.loc, end = ins.loc + ins.length, type = "insertion")
    } else {
      all_mutations <- strsplit(x = mutation, split = "_")[[1]]
      locations <- as.numeric(sapply(X = all_mutations, FUN = gsub, pattern = "([A-Z,a-z]+)", replacement = ""))
      locations <- locations[!is.na(locations)]
      locations <- data.frame(start = locations + segment.df[segment.df$segment == segment, "start"],
                              type = "nonsynonymous")
      if (any(all_mutations == "synonymous")) {
        locations[which(all_mutations == "synonymous") - 1, "type"] <- "synonymous"
      }
      if (any(all_mutations == "noncoding")) {
        locations[which(all_mutations == "noncoding") - 1, "type"] <- "noncoding"
      }
      locations$end <- locations$start
      mutation.df[[length(x = mutation.df) + 1]] <- locations
    }
  }
  return.df <- do.call(what = rbind, args = c(mutation.df, deletion.df, insertion.df))
  return.df$cell <- cell
  return(return.df)
}

# Create long data.frame of segment presence for each cell
#
# @param data data.frame containing segment presence information
# @param segment.df data.frame containing the coordinates of each segment
# @param pattern prefix pattern to match segment presence columns
#
# @return Returns a data.frame where each row contains a detected segment, its
#         positional information, and cell index

ExtractSegmentPresence <- function(data, segment.df, pattern) {
  segment.df.per.cell <- lapply(X = 1:nrow(x = data), FUN = function(x) {
    segment.df$cell <- x
    segment.df$present <- as.logical(x = unlist(x = data[x, paste0(pattern, segment.df$segment)]))
    segment.df <- segment.df[segment.df$present, ]
  })
  segment.df.per.cell <- do.call(what = rbind, args = segment.df.per.cell)
  segment.df.per.cell$present <- NULL
  return(segment.df.per.cell)
}

# Main ggplot2 theme definition
PacBioTheme <- function () {
  pacbio.theme <- theme(
    axis.ticks.y=element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size = 12),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(),
    plot.margin = margin(0, 20, 90, 0),
    legend.key = element_rect(fill = NA),
    legend.title=element_text(size = 15),
    legend.text=element_text(size = 12)
  )
  return(pacbio.theme)
}

### Main plotting function #####################################################

PacBioPlot <- function(
    data,
    segment.df,
    color.segments.by = "percent_viral_UMIs",
    segment.color = '#0072B2',
    order.segments.by = "percent_viral_UMIs",
    order.name = "% mRNA from flu",
    box = 'percent_supernatant',
    box.name ='% supernatant',
    box.color.high = '009E73',
    ncol = 1,
    arrow_height = 0.325, # vertical space per arrow in inches
    arrow_frac_height = 0.55 # arrow takes up this much of available height
) {
  ncells <- nrow(x = data)
  segments <- segment.df$segment
  segment.df.cells <- ExtractSegmentPresence(data = data, segment.df = segment.df, pattern = "present_")
  mutations.df <- ParseMutations(data = data, segments = segments, segment.df = segment.df, pattern = "mutations_")
  segment.df.cells <- segment.df.cells[segment.df.cells$cell <= ncells, ]
  mutations.df <- mutations.df[mutations.df$cell <= ncells, ]
  segment.df.cells[['segment_coloring']] <- dat[[color.segments.by]][segment.df.cells$cell]
  if (!is.null(x = box)) {
    if (length(x = box) != length(x = box.name)) {
      stop("Please provide equal length vectors for the parameters box and box.name.")
    }
    box.col.names <- paste0('box_', 1:length(x = box))
    for (i in 1:length(x = box)) {
      segment.df.cells[[box.col.names[i]]] <- dat[[box[i]]][segment.df.cells$cell]
    }
  }
  # reorder cells
  segment.df.cells[['segment_ordering']] <- dat[[order.segments.by]][segment.df.cells$cell]
  cell.remap <- 1:ncells
  names(x = cell.remap) <- order(unique(x = segment.df.cells[, c("cell", "segment_ordering")])$segment_ordering, decreasing = TRUE)
  segment.df.cells$cell_ordered <- sapply(X = segment.df.cells$cell, FUN = function(x) cell.remap[as.character(x = x)])
  mutations.df$cell_ordered <- sapply(X = mutations.df$cell, FUN = function(x) cell.remap[as.character(x = x)])
  # dummy data.frame for plotting missing segments
  missing.segment.df <- data.frame(y = 1:ncells)
  indel.df <- mutations.df[mutations.df$type %in% c("insertion", "deletion"), ]
  indel.df$type <- droplevels(indel.df$type)
  mutations.df <- mutations.df[!mutations.df$type %in% c("insertion", "deletion"), ]
  mutations.df$type <- droplevels(x = mutations.df$type)
  box.labels <- unique(x = segment.df.cells[, c(box.col.names, "cell_ordered")])
  box.labels[, 1:length(x = box.col.names)] <- apply(box.labels[, 1:length(x = box.col.names), drop = FALSE], MARGIN = 2, FUN = round, digits = 1)
  cell.cols <- split(1:ncells, ceiling(seq_along(1:ncells)/(ncells/ncol)))
  plots <- lapply(X = cell.cols, FUN = function(x) {
    ncells.plot <- length(x = x)
    segment.df.cells.sub <- segment.df.cells[segment.df.cells$cell_ordered %in% x, ]
    segment.df.cells.sub$cell_ordered <- abs(segment.df.cells.sub[segment.df.cells.sub$cell_ordered %in% x, "cell_ordered"] - (max(x) + 1))
    missing.segment.df.sub <- missing.segment.df[missing.segment.df$y %in% x, ,drop=FALSE]
    missing.segment.df.sub$y <- abs(missing.segment.df.sub$y - (max(x) + 1))
    mutations.df.sub <- mutations.df[mutations.df$cell_ordered %in% x, ]
    mutations.df.sub$cell_ordered <- abs(mutations.df.sub[mutations.df.sub$cell_ordered %in% x, "cell_ordered"] - (max(x) + 1))
    indel.df.sub <- indel.df[indel.df$cell_ordered %in% x, ]
    indel.df.sub$cell_ordered  <- abs(indel.df.sub[indel.df.sub$cell_ordered %in% x, "cell_ordered"] - (max(x) + 1))
    box.labels.sub <- box.labels[box.labels$cell_ordered %in% x, ]
    box.labels.sub$cell_ordered <- abs(box.labels.sub[box.labels.sub$cell_ordered %in% x, "cell_ordered"] - (max(x) + 1))
    plot <- ggplot(segment.df.cells.sub) +
      geom_segment(data = missing.segment.df.sub, aes_string(y = "y", yend = "y", x = 0, xend= max(segment.df$end)), color = "grey") +
      geom_gene_arrow(aes_string(xmin = "start", xmax = "end", y = "cell_ordered"),
                      arrow_body_height = unit(arrow_frac_height * arrow_height, 'in'),
                      arrowhead_height = unit(arrow_frac_height * arrow_height, 'in'),
                      arrowhead_width = unit(0.3 * arrow_height, 'in'),
                      fill = segment.color,
                      size = 0.4) +
      geom_segment(data = indel.df.sub, aes_string(x = "start", xend = "end", y = "cell_ordered", yend = "cell_ordered", color = "type"), size = 2.5) +
      scale_color_manual(values = c('#CC79A7', '#F0E442'), name = "Indel Class", drop = FALSE, guide = guide_legend(order = 2)) +
      new_scale_color() +
      geom_point(data = mutations.df.sub, aes_string(x = "start", y = "cell_ordered", fill = "type"), color = "black", size = 3.5, pch = 21) +
      scale_fill_manual(values = c("#999999", '#E69F00', "#009E73"), name = "Mutation Class", drop = FALSE, guide = guide_legend(order = 1)) +
      new_scale_fill()
    x.pos <- NULL
    if (!is.null(x = box)) {
      for (i in 1:length(x = box)) {
        x.pos[i] <- -500 - ((i-1) * 500)
        box.labels.sub[is.na(box.labels.sub)] <- "NA"
        plot <- plot + geom_tile(aes_string(x = x.pos[i], y = "cell_ordered", width = 400, height = arrow_frac_height + 0.1, fill = box.col.names[i]), color = "black") +
          geom_text(data = box.labels.sub, aes_string(x = x.pos[i], y = "cell_ordered", label = box.col.names[i]), size = 3) +
          scale_fill_gradient(low = 'white', high = box.color.high[i], name = box.name[i], limits = c(min(segment.df.cells[[box.col.names[i]]]), max(segment.df.cells[[box.col.names[i]]]))) +
          new_scale_fill()
        if (length(x = box) > 1) {
          plot <- plot + coord_cartesian(clip = 'off') +
            annotate(geom = "text", x = x.pos[i], y = 0, label = box.name[i], angle = 90, size = 4, hjust = 1)
        }
      }
    }
    if (length(x = box) > 1) {
      plot <- plot +
        scale_x_continuous("", breaks = c(x.pos, segment.df$start + segment.df$len/2), labels = c(rep("", times = length(x = box)), gsub(pattern = "flu", replacement = "", x = segment.df$segment)), limits = c(min(x.pos)-200, max(segment.df$end)), expand = c(0.01, 0.01)) +
        scale_y_continuous("", breaks = 1:ncells.plot, labels = paste0('cell ', rev(x))) +
        coord_cartesian(clip = "off", ylim = c(0.3, ncells.plot + 0.5), expand = FALSE) +
        PacBioTheme()
    } else {
      plot <- plot +
        scale_x_continuous("", breaks = c(x.pos, segment.df$start + segment.df$len/2), labels = c(box.name, gsub(pattern = "flu", replacement = "", x = segment.df$segment)), limits = c(min(x.pos)-200, max(segment.df$end)), expand = c(0.01, 0.01)) +
        scale_y_continuous("", breaks = 1:ncells.plot, labels = paste0('cell ', rev(x)), expand = c(0.01, 0.01), limits = c(0.5, ncells.plot + 0.5)) +
        PacBioTheme() + theme(plot.margin = margin(0, 20, 0, 0))
    }
  })
  wrap_plots(plots) + plot_layout(guides = 'collect')
}

### Setup ######################################################################
# Read in the data and generate influenza segment data.frame with positioning 
# information.

dat <- read.csv(file = csv.filename)
dat$percent_viral_UMIs <- dat$frac_viral_UMIs * 100
dat$percent_supernatant <- dat$freq_supernatant * 100

segment.df <- data.frame(
  segment = c("fluPB2", "fluPB1", "fluPA", "fluHA", "fluNP", "fluNA", "fluM", "fluNS"),
  len = c(2341, 2341, 2233, 2035, 1565, 1735, 1027, 890)
)
segment.df$end <- cumsum(segment.df$len)
segment.df$start <- segment.df$end - segment.df$len + 1

### Generate/Save plot #########################################################
plot <- PacBioPlot(
  dat,
  segment.df,
  box = c("percent_viral_UMIs", "percent_supernatant"),
  box.name = c("% mRNA from flu", "% supernatant"),
  box.color.high = c("#009E73", "#E69F00"),
  order.segments.by = "percent_viral_UMIs",
  order.name = "% supernatant",
  color.segments.by = 'percent_supernatant',
  ncol = 1
)
ggsave(plot, filename = plot.filename, width = 14, height = 30)
