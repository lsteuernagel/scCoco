#' retrieve aba ids for a vector of genes
#'
#' retreives all ids associated with provided ids. Genes that are not included in aba will be silently dropped!
#'
#' @param genes a vector of gene symbols
#' @param aba_gene_to_id a dataframe ('id','gene_symbol'), a filepath to a valid dataframe or NULL. if NULL: will use cocoframer to query
#' @return  vector with ids and genes as names
#'
#' @export
#'
#' @import cocoframer
#'
#' @importFrom data.table fwrite
#'
#

aba_ids =  function(genes, aba_gene_to_id =NULL){

  ### get api ids

  if(!is.null(aba_gene_to_id)){
    # query from provided or loaded file
    if(is.character(aba_gene_to_id)){
      if(file.exists(aba_gene_to_id)){
        aba_gene_to_id = data.table::fread(aba_gene_to_id,data.table = FALSE)
      }else{
        stop("Detected string to cached file provided via aba_gene_to_id but cannot find file. Please provide a correct path, a valid data.frame or set to NULL to query via the API.")
      }
    }else if(is.data.frame(aba_gene_to_id)){
      if(colnames(aba_gene_to_id)[1] != "id" | colnames(aba_gene_to_id)[2] != "gene_symbol"){
        stop("Detected string to cached file provided via aba_gene_to_id but cannot find file. Please provide a correct path, a valid data.frame or set to NULL to query via the API.")
      }
    }else{
      stop("Not a dataframe or string.  Please provide a correct path, a valid data.frame or set to NULL to query via the API.")
    }

    #  get ids
    all_ids = as.character(aba_gene_to_id$id[aba_gene_to_id$gene_symbol %in% genes])
    names(all_ids) = aba_gene_to_id$gene_symbol[aba_gene_to_id$gene_symbol %in% genes]

  }else{
    # api query
    all_ids = character(0)
    for(i in 1:length(genes)){
      current_id = cocoframer::get_gene_aba_ish_ids(genes[i])
      names(current_id) = rep(genes[i],length(current_id))
      all_ids = c(all_ids,current_id)
    }
  }
  return(all_ids)
}

#' Retrieve aba expression data per voxel for given ids
#'
#' Retrieves expression data per voxel.
#'
#' @param ids a vector of aba ids
#' @param return should the function return list or matrix (see below)
#' @param aba_ish_list a list of 3d arrays or a string to a cached version (.rds)
#' @param aba_ish_matrix a voxel x gene matrix or a string to a cached version (table)
#' @return  a list or matrix depending on return
#'
#' @export
#'
#' @import cocoframer
#'
#' @importFrom data.table fwrite
#'

coco_expression = function(ids,return="matrix",aba_ish_list=NULL,aba_ish_matrix = NULL){
  stop_convert=FALSE
  # if nothing is provided query API:
  if(is.null(aba_ish_list) & is.null(aba_ish_matrix)){
    aba_ish_res =list()
    message("Querying ABA API with ids. May take some time.")
    for(i in 1:length(ids)){
      current_id = ids[i]
      temp <- tryCatch({
        cocoframer::get_aba_ish_data(api_id=current_id, values = "energy")
      },
      error=function(cond) {
        message("Query failed for:",current_id, " with error: ",cond)
        return(NULL)
      })
      aba_ish_res[[current_id]] = temp
    }
    # else load or use matrix if provided
  }else if(is.null(aba_ish_list) & !is.null(aba_ish_matrix)){
    # load matrix
    if(is.character(aba_ish_matrix)){
      if(file.exists(aba_ish_matrix)){
        aba_ish_res = data.table::fread(aba_ish_matrix,data.table = FALSE,header = TRUE)
      }else{
        stop("Detected string to cached file provided via aba_ish_matrix but cannot find file. Please provide a correct path, a valid data.frame or set to NULL to query via the API.")
      }
    }else if(is.matrix(aba_ish_matrix) | is.data.frame(aba_ish_matrix)){
      aba_ish_res = aba_ish_matrix
      stop_convert =TRUE
    }else{
      stop("Not a matrix or string.  Please provide a correct path, a valid matrix or set to NULL to query via the API.")
    }
    aba_ish_res = aba_ish_res[,ids,drop=FALSE]
    # else load or use rds list if provided
  }else if(!is.null(aba_ish_list) & is.null(aba_ish_matrix)){
    # load list
    if(is.character(aba_ish_list)){
      if(file.exists(aba_ish_list)){
        aba_ish_res = readRDS(aba_ish_list)
      }else{
        stop("Detected string to cached file provided via aba_ish_list but cannot find file. Please provide a correct path, a valid data.frame or set to NULL to query via the API.")
      }
      aba_ish_res = aba_ish_list
    }else{
      stop("Not a matrix or string.  Please provide a correct path, a valid matrix or set to NULL to query via the API.")
    }
    if(length(ids)==1){
      aba_ish_res=list(aba_ish_res[ids])
      names(aba_ish_res) = ids
    }else{
      aba_ish_res = aba_ish_res[ids]
    }
  }
  if(return == "matrix" & !stop_convert){
    if(length(aba_ish_res)>1){
      aba_ish_res = as.data.frame(do.call(cbind,aba_ish_res))
    }else{
      aba_ish_res = as.data.frame(aba_ish_res)
    }
  }
  return(aba_ish_res)
}
