% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readCOM.R
\name{readCOM}
\alias{readCOM}
\title{Loads multiple or a single Center Of Mass (COM) or Centroid .xyz file(s) and 
Converts to 2D x and y coordinates based on a viewing matrix}
\usage{
readCOM(directory, pmat, seperator = " ", extension = ".xyz")
}
\arguments{
\item{directory}{Path containing the input .xyz file(s)}

\item{pmat}{A 4 x 4 Viewing Matrix used to transform the 3D coordinates to 2D}

\item{seperator}{how the columns are seperated. Example: " " for space, "," for csv files, "\t" for tab-delim files.}
}
\value{
list with transformed coordinates and min/max + length of files loaded info
}
\description{
NB*** Ensure the extension of the files is .xyz, contains no header and is column ordered X, Y, Z.
}
