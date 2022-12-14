% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression_functions.R
\name{aba_ids}
\alias{aba_ids}
\title{Retrieve aba ids for a vector of genes}
\usage{
aba_ids(genes, aba_gene_to_id = NULL)
}
\arguments{
\item{genes}{a vector of gene symbols}

\item{aba_gene_to_id}{a dataframe ('id','gene_symbol'), a filepath to a valid dataframe that can be read with data.table::fread or NULL. if NULL: will use cocoframer::get_gene_aba_ish_ids to query}
}
\value{
vector with (mba) ids and genes as names
}
\description{
Retrieves all ids associated with provided ids. Genes that are not included in aba will be silently dropped!
}
