#' Median rank across a matrix of gene expression
#'
#'
#' @param expression_matrix expression_matrix a matrix or dataframe with rows as voxels and columns as genes
#' @param subset numeric vector of ids to subset matrix to. defaults to NULL (no subset)
#' @param convert do 1 - (expression_matrix_rank_median / nrow(expression_matrix)) ?
#' @return  vector with median ranks per voxel and numeric ids as names
#'
#' @export
#'
#'
#
median_rank = function(expression_matrix,subset=NULL,convert =TRUE){
  if(!is.null(subset)){
    expression_matrix = expression_matrix[subset,]
  }
  # make ranks
  expression_matrix_rank = nrow(expression_matrix) - apply(expression_matrix,2,rank,ties="min")
  expression_matrix_rank_median = apply(expression_matrix_rank,1,median)
  names(expression_matrix_rank_median) = 1:nrow(expression_matrix_rank)
  # convert
  expression_matrix_rank_median = 1 - (expression_matrix_rank_median / nrow(expression_matrix))
  #return(expression_matrix_rank_median)
  return(expression_matrix_rank_median)
}

#' Annotate regions via voxels using ABA ISH expression data and an enrichment function.
#'
#' Function uses column names of aba ish objects as provided when queried using cocoframer.
#'
#' @param expression_matrix e
#' @param target_structure_id defaults to 1097 for hypothalamus. set to 997 for root!
#' @param target_level
#' @param enrich_function
#' @param aba_ccf_grid_annotation
#' @param  mba_ontology_flatten
#' @return  vector with median ranks per voxel and numeric ids as names
#'
#' @export
#'
#' @import dplyr data.table
#'
#'

region_annotation = function(expression_matrix,target_structure_id = "1097",target_level = "8",enrich_function = "median_rank", aba_ccf_grid_annotation = NULL ,mba_ontology_flatten= NULL){

  # get ccf grid as 3d array
  if(!is.null(aba_ccf_grid_annotation)){
    # query from provided or loaded file
    if(is.character(aba_ccf_grid_annotation)){
      if(file.exists(aba_ccf_grid_annotation)){
        aba_ccf_grid_annotation = readRDS(aba_ccf_grid_annotation)
      }else{
        stop("Detected string to cached file provided via aba_ccf_grid_annotation but cannot find file. Please provide a correct path, a valid array or set to NULL to query via the API.")
      }
    }else if(is.array(aba_gene_to_id)){

    }else{
      stop("Not a 3d array or string.  Please provide a correct path, a valid array or set to NULL to query via the API.")
    }
  }else{
    # api query
    aba_ccf_grid_annotation = cocoframer::get_ccf_grid_annotation()
  }

  # get region info
  if(!is.null(mba_ontology_flatten)){
    # query from provided or loaded file
    if(is.character(mba_ontology_flatten)){
      if(file.exists(mba_ontology_flatten)){
        mba_ontology_flatten = data.table::fread(mba_ontology_flatten,data.table = FALSE,header = TRUE)
      }else{
        stop("Detected string to cached file provided via mba_ontology_flatten but cannot find file. Please provide a correct path, a valid data.frame or set to NULL to query via the API.")
      }
    }else if(is.data.frame(mba_ontology_flatten)){

    }else{
      stop("Not a data.frame or string.  Please provide a correct path, a valid array or set to NULL to query via the API.")
    }
  }else{
    # api query
    mba_ontology = cocoframer::get_mba_ontology()
    mba_ontology_flatten = cocoframer::flatten_mba_ontology(mba_ontology)
  }
  # set some ids to character
  mba_ontology_flatten$parent_structure_id = as.character(mba_ontology_flatten$parent_structure_id)
  mba_ontology_flatten$atlas_id = as.character(mba_ontology_flatten$atlas_id)
  mba_ontology_flatten$id = as.character(mba_ontology_flatten$id)

  # find all children for selected node to subset
  # hypothalamus specific structures:
  mba_ontology_flatten_edges = data.frame(from=as.character(mba_ontology_flatten$parent_structure_id), to = as.character(mba_ontology_flatten$id),stringsAsFactors = F)
  target_structure_children = find_children(mba_ontology_flatten$id[mba_ontology_flatten$id %in% target_structure_id],edges = mba_ontology_flatten_edges)

  # take voxel anno and mba_flattened
  full_annotation = data.frame(voxel_id = 1:nrow(expression_matrix),structure_id = as.character(aba_ccf_grid_annotation))
  full_annotation = dplyr::left_join(full_annotation,mba_ontology_flatten,by=c("structure_id"="id"))
  # which voxels to keep
  expression_matrix_subset = expression_matrix[which(full_annotation$structure_id %in% target_structure_children),]

  # run scoring function:
  enrich_function = match.fun(enrich_function)
  score_subset = enrich_function(expression_matrix_subset,convert=TRUE)

  # result
  all_target_voxels = full_annotation[which(full_annotation$structure_id %in% target_structure_children),]
  all_target_voxels$score = score_subset

  # convert to matrix per region
  all_target_voxels_summarize = all_target_voxels %>% dplyr::group_by(structure_id,name) %>% dplyr::slice_max(order_by = score,n=topn) %>%
    dplyr::summarise(region_median = median(score)) %>% suppressMessages() %>% as.data.frame()

  # further summarize to specific level
  mba_ontology_flatten$st_level = as.numeric(mba_ontology_flatten$st_level)
  structures_on_target_level = find_children(nodes=target_structure_id,edges=mba_ontology_flatten_edges[mba_ontology_flatten$st_level <= as.numeric(target_level),])
  structures_on_target_level = structures_on_target_level[structures_on_target_level %in% mba_ontology_flatten$id[mba_ontology_flatten$st_level == as.numeric(target_level)]]
  structures_on_target_level_withchildren = sapply(structures_on_target_level,find_children,edges = mba_ontology_flatten_edges)
  structures_on_target_level_withchildren = sapply(names(structures_on_target_level_withchildren),function(name,list){c(name,list[[name]])},list=structures_on_target_level_withchildren)
  structures_on_target_level_withchildren_df = lapply(names(structures_on_target_level_withchildren),function(name,list){df = data.frame(subnodes = list[[name]], topnode = name);return(df)},list=structures_on_target_level_withchildren)
  structures_on_target_level_withchildren_df = do.call(rbind.data.frame,structures_on_target_level_withchildren_df) %>% dplyr::filter(topnode %in% structures_on_target_level)
  ## add to voxel df
  all_target_voxels =  dplyr::left_join(all_target_voxels,structures_on_target_level_withchildren_df,by=c("structure_id"="subnodes"))
  all_target_voxels$topnode[is.na(all_target_voxels$topnode)] = all_target_voxels$structure_id[is.na(all_target_voxels$topnode)]
  all_target_voxels = dplyr::left_join(all_target_voxels,mba_ontology_flatten %>% dplyr::select(id, topname = name),by=c("topnode"="id"))
  # further sum
  all_target_voxels_summarize_targetlevel = all_target_voxels %>% dplyr::group_by(topnode,topname) %>% dplyr::slice_max(order_by = score,n=topn) %>%
    dplyr::summarise(region_median = median(score)) %>% suppressMessages() %>% as.data.frame()

  # make annotated version of expression matrix
  expression_matrix_subset_res = cbind(voxel_id = rownames(expression_matrix_subset),expression_matrix_subset)

  # group for output object
  return_list = list(
   # return matrix with scores
    genes_per_voxel = expression_matrix_subset_res,
   # return all_target_voxels
    scores_per_voxel_annotated = all_target_voxels %>% as.data.frame(),
   # return summary per ..
   scores_per_leaf_region = all_target_voxels_summarize,
   # return per target region
   scores_per_target_level_region = all_target_voxels_summarize_targetlevel
  )
  return(return_list)
}

#' Recursive helper function for graph traversal
#'
#' @param nodes starting node(s)
#' @param edges edgelist (data.frame: from and to)
#' @return  vector with all_children from nodes
#'
#' @export
#'
#'

find_children = function(nodes,edges){
  current_children = edges$to[edges$from %in% nodes]
  #print(paste0(current_children,collapse = "|"))
  if(length(current_children)>0){
    all_children = c(current_children,find_children(current_children,edges))
  }else{
    all_children = current_children
  }
  return(all_children)
}

## summarise to dataframe function:
  ## need function that returns likely and other clusters as data.frame

## need wrapper functio(s)
  # execute above pipeline for sets of marker genes and store (intermediate) results
  # need to allow confidence measure based on detected genes




