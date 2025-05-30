% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{venn}
\alias{venn}
\title{Create Venn Diagrams with Customizable Appearance}
\usage{
venn(
  ...,
  out_name = NULL,
  PPT = NULL,
  label = NULL,
  label.cex = 8,
  text.cex = 8,
  fill.alpha = 1,
  fill.col = c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0"),
  stroke.alpha = 1,
  stroke.size = 1,
  mar = 1,
  show.elements = FALSE,
  label.sep = "\\n",
  auto.scale = FALSE,
  x.lim = NULL,
  y.lim = NULL,
  add_border = FALSE,
  hjust = NULL,
  stroke.linetype = "longdash",
  col.stroke = "black",
  show.percentage = FALSE,
  w = 7,
  h = 5
)
}
\arguments{
\item{...}{One or more vectors containing elements to compare (max 4)}

\item{out_name}{Character, output filename prefix (without extension), NULL means don't save}

\item{PPT}{Character, if not NULL, saves plot as PowerPoint with this name}

\item{label}{Character vector, custom names for each input set}

\item{label.cex}{Numeric, size of set labels, default is 8}

\item{text.cex}{Numeric, size of intersection counts, default is 8}

\item{fill.alpha}{Numeric, transparency of fill colors (0-1), default is 1}

\item{fill.col}{Character vector, fill colors for each set}

\item{stroke.alpha}{Numeric, transparency of stroke lines (0-1), default is 1}

\item{stroke.size}{Numeric, thickness of stroke lines, default is 1}

\item{mar}{Numeric, margin around plot, default is 1}

\item{show.elements}{Logical, whether to show elements instead of counts, default is FALSE}

\item{label.sep}{Character, separator for element labels when show.elements=TRUE, default is newline}

\item{auto.scale}{Logical, whether to auto-scale diagram, default is FALSE}

\item{x.lim}{Numeric vector, x-axis limits, NULL means automatic}

\item{y.lim}{Numeric vector, y-axis limits, NULL means automatic}

\item{add_border}{Logical, whether to add plot border, default is FALSE}

\item{hjust}{Numeric, horizontal justification of labels, NULL means automatic}

\item{stroke.linetype}{Character, line type for set boundaries, default is "longdash"}

\item{col.stroke}{Character, color for stroke lines, default is "black"}

\item{show.percentage}{Logical, whether to show percentages instead of counts, default is FALSE}

\item{w}{Numeric, output width in inches, default is 7}

\item{h}{Numeric, output height in inches, default is 5}
}
\value{
Returns a ggplot object containing the Venn diagram. For >4 sets, returns intersection matrix.
}
\description{
Generates Venn diagrams to visualize overlaps between multiple sets of elements.
Supports 2-4 sets and provides extensive customization options for appearance.
}
\examples{
set1 <- letters[1:10]
set2 <- letters[5:15]
venn(set1, set2, label = c("Group A", "Group B"))

}
