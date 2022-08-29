#' Median rank across a matrix of gene expression
#'
#'
#' @param expression_matrix expression_matrix a matrix or dataframe with rows as voxels and columns as genes
#' @param subset numeric vector of ids to subset matrix to. defaults to NULL (no subset)
#' @param convert do 1 - (expression_matrix_rank_median / nrow(expression_matrix)) ?
#' @param weights weights for ids
#' @return  vector with median ranks per voxel and numeric ids as names
#'
#' @importFrom stats median
#'
#' @export
#'
#'
#
median_rank = function(expression_matrix,subset=NULL,convert =TRUE,weights=NULL){
  if(!is.null(subset)){
    expression_matrix = expression_matrix[subset,]
  }
  # make ranks
  expression_matrix_rank = nrow(expression_matrix) - apply(expression_matrix,2,rank,ties="min")
  if(!is.null(weights)){
    weights = ceiling(weights)
    repeated_colnames = rep(colnames(expression_matrix_rank),weights)
    expression_matrix_rank = expression_matrix_rank[,repeated_colnames]
  }
  expression_matrix_rank_median = apply(expression_matrix_rank,1,stats::median)
  names(expression_matrix_rank_median) = 1:nrow(expression_matrix_rank)
  # convert
  if(convert){
    expression_matrix_rank_median = 1 - (expression_matrix_rank_median / nrow(expression_matrix))
  }
  #return(expression_matrix_rank_median)
  return(expression_matrix_rank_median)
}

#' Annotate regions via voxels using ABA ISH expression data and an enrichment function.
#'
#' Function uses column names of aba ish objects as provided when queried using cocoframer.
#'
#' @param expression_matrix expression_matrix a matrix or dataframe with rows as voxels and columns as genes
#' @param weights weights for ids (columns in expression_matrix). provide null for unweighted
#' @param target_structure_id defaults to 997 for hypothalamus. set to 1097 for hypothalamus.!
#' @param exclude_ids mba structure ids to exclude. defaults to subfornical organ
#' @param target_level number st_level in mba data (e.g. level 8 to provide a summary on this level)
#' @param enrich_function a function (defaults to 'median_rank') for enrichment
#' @param topn topn voxels to annotate region. defaults to 10. set to inf if not filtering for averaging is desired
#' @param aba_ccf_grid_annotation a 3d array from cocoframer::get_ccf_grid_annotation() or a filepath to a rds containing a cached version. NULL will query the API
#' @param  mba_ontology_flatten a dataframe from  cocoframer::get_mba_ontology() + cocoframer::flatten_mba_ontology() or a filepath to a tsv containing a cached version. NULL will query the API
#' @return  vector with median ranks per voxel and numeric ids as names
#'
#' @export
#'
#' @import dplyr magrittr
#'
#' @importFrom data.table fwrite
#'
#'

region_annotation = function(expression_matrix,weights=NULL,target_structure_id = "997",exclude_ids = c("338"),target_level = "8",enrich_function = "median_rank",topn=10, aba_ccf_grid_annotation = NULL ,mba_ontology_flatten= NULL){

  # get ccf grid as 3d array
  if(!is.null(aba_ccf_grid_annotation)){
    # query from provided or loaded file
    if(is.character(aba_ccf_grid_annotation)){
      if(file.exists(aba_ccf_grid_annotation)){
        aba_ccf_grid_annotation = readRDS(aba_ccf_grid_annotation)
      }else{
        stop("Detected string to cached file provided via aba_ccf_grid_annotation but cannot find file. Please provide a correct path, a valid array or set to NULL to query via the API.")
      }
    }else if(is.array(aba_ccf_grid_annotation)){
      # do nothing
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

  # exclude
  target_structure_children = target_structure_children[! as.character(target_structure_children) %in% as.character(exclude_ids)]

  # take voxel anno and mba_flattened
  full_annotation = data.frame(voxel_id = 1:nrow(expression_matrix),structure_id = as.character(aba_ccf_grid_annotation))
  full_annotation = dplyr::left_join(full_annotation,mba_ontology_flatten,by=c("structure_id"="id"))
  # which voxels to keep
  expression_matrix_subset = expression_matrix[which(full_annotation$structure_id %in% target_structure_children),]

  # run scoring function:
  enrich_function = match.fun(enrich_function)
  score_subset = enrich_function(expression_matrix_subset,convert=TRUE,weights=weights)

  # result
  all_target_voxels = full_annotation[which(full_annotation$structure_id %in% target_structure_children),]
  all_target_voxels$score = score_subset

  # convert to matrix per region
  all_target_voxels_summarize = all_target_voxels %>% dplyr::group_by(structure_id,name) %>% dplyr::slice_max(order_by = score,n=topn) %>%
    dplyr::summarise(region_median = median(score)) %>% suppressMessages() %>% as.data.frame()

  # further summarize to specific level
  mba_ontology_flatten$st_level = as.numeric(mba_ontology_flatten$st_level)

  # this version also includes cases where the tree jumps over the target lvele , i.e. a parent only having children on alevel below (the these will be included too!)
  structures_selected = find_children(nodes=target_structure_id,edges=mba_ontology_flatten_edges)
  structures_above_target_level = structures_selected[structures_selected %in% mba_ontology_flatten$id[mba_ontology_flatten$st_level < as.numeric(target_level)]]
  children_structures_above_target_level =  mba_ontology_flatten_edges$to[mba_ontology_flatten_edges$from %in%structures_above_target_level]
  structures_on_target_level = children_structures_above_target_level[children_structures_above_target_level %in% mba_ontology_flatten$id[mba_ontology_flatten$st_level >= as.numeric(target_level)]]

  # get children of these into a list
  structures_on_target_level = structures_on_target_level[structures_on_target_level %in% mba_ontology_flatten$id[mba_ontology_flatten$st_level >= as.numeric(target_level)]]
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

#' Annotate sets of genes (cluster Markers) with mouse brain regions
#'
#' TODO: need to explain all results in details
#'
#' @param gene_set gene set(s). A list of character vectors
#' @param min_ids for quality control (pct_of_genes): How many ids should be included in the calculation. else pct_of_genes will be 0
#' @inheritParams coco_expression
#' @inheritParams aba_ids
#' @param max_ids_to_include maximum of ids that will bes used fr enrichment. defaults to inf
#' @param ... parameters passed down to \code{\link{region_annotation}}
#' @return  a list with results
#'
#' @export
#'
#'

findRegions_genesets = function(gene_set,gene_set_weights=NULL,min_ids = 0, aba_gene_to_id =NULL, aba_ish_matrix =NULL,max_ids_to_include = Inf,...){

  # check whether gene_set is generally valid input
  if(is.list(gene_set)){
    if(any(!sapply(gene_set,is.character))){
      stop("Please provide a character vector or a list of character vectors with mouse gene symbols" )
    }
  }else{
    if(is.character(gene_set)){
      gene_set = list(gene_set = gene_set)
    }else{
      stop("Please provide a character vector or a list of character vectors with mouse gene symbols" )
    }
  }

  if(!is.null(gene_set_weights)){
    if(any(!sapply(gene_set_weights,is.numeric))){
      stop("Please provide a numeric vector or a list of numeric vectors with weights for each gene set" )
    }
  }

  ## init lists for results:
  genes_per_voxel_list =list()
  pct_of_genes_list =list()
  scores_per_voxel_annotated_list =list()
  scores_per_leaf_region_list =list()
  scores_per_target_level_region_list =list()
  queried_genes_list = list()
  ids_length_list = list()
  pct_of_weights_list =list()

  ## run over all entries
  for(i in 1:length(gene_set)){

    # get
    ids_to_query = aba_ids(gene_set[[i]],aba_gene_to_id = aba_gene_to_id)
    queried_genes_list[[names(gene_set)[i]]] = unique(names(ids_to_query))
    ids_length = length(unique(names(ids_to_query)))
    ids_length_list[[names(gene_set)[i]]] = ids_length
    if(ids_length < min_ids){
      ids_length = 0
    }
    if(ids_length > max_ids_to_include){
      ids_to_query = ids_to_query[1:max_ids_to_include]
    }
    pct_of_genes_list[[names(gene_set)[i]]] = ids_length / length(unique(gene_set[[i]]))
    if(length(ids_to_query) < 2){
      message("Cannot find any (or only 1) ids for gene in set ",names(gene_set)[i])
      pct_of_weights_list[[names(gene_set)[i]]] = NA
    }else{
      message("Running ",names(gene_set)[i]," with ",length(ids_to_query)," ids to query.")
      # get expression
      aba_expression = coco_expression(ids = ids_to_query,aba_ish_matrix = aba_ish_matrix,return = "matrix")
      # optionally apply weights
      if(!is.null(gene_set_weights)){
        weights = gene_set_weights[[i]]
        weights_per_id = weights[match(names(ids_to_query),gene_set[[i]])]
        pct_of_weights_list[[names(gene_set)[i]]] = sum(weights_per_id) / sum(weights)
        if(length(weights_per_id) == ncol(aba_expression)){
          temp_res = region_annotation(expression_matrix = aba_expression, weights = weights_per_id,...)
        }else{
          pct_of_weights_list[[names(gene_set)[i]]] = 1
          temp_res = region_annotation(expression_matrix = aba_expression, ...)
          message("Warning: Cannot apply weights: different lengths.")
        }
      }else{
        pct_of_weights_list[[names(gene_set)[i]]] = 1
        temp_res = region_annotation(expression_matrix = aba_expression, ...)
      }
      # save in individual lists
      genes_per_voxel_list[[names(gene_set)[i]]] = temp_res$genes_per_voxel
      scores_per_voxel_annotated_list[[names(gene_set)[i]]] = temp_res$scores_per_voxel_annotated
      scores_per_leaf_region_list[[names(gene_set)[i]]] = temp_res$scores_per_leaf_region
      scores_per_target_level_region_list[[names(gene_set)[i]]] = temp_res$scores_per_target_level_region
    }
  }

  # bind matrices
  templist = lapply(scores_per_leaf_region_list,function(x){x[,3]})
  tempmat = do.call(cbind,templist)
  scores_per_leaf_region_all = cbind(scores_per_leaf_region_list[[1]][,1:2],tempmat)
  templist = lapply(scores_per_target_level_region_list,function(x){x[,3]})
  tempmat = do.call(cbind,templist)
  scores_per_target_level_region_all = cbind(scores_per_target_level_region_list[[1]][,1:2],tempmat)

  #### QC
  # pct of genes as QC metric:
  pct_of_genes_all = as.data.frame(do.call(rbind,pct_of_genes_list))
  colnames(pct_of_genes_all) = "pct_of_genes"
  pct_of_genes_all$number_ids = do.call(rbind,ids_length_list)[,1]
  # overlap of lists
  intersect_matrix = intersect_from_list(gene_set)
  geneset_binary_dist = as.matrix(dist(t(intersect_matrix),method = "binary")) # return this!
  # get the smallest dist
  shortest_dists = apply(geneset_binary_dist,1,function(x){
    second_lowest_dist = sort(x)[2]
    names(second_lowest_dist) = names(x[which(x==second_lowest_dist)])[1]
    second_lowest_dist
  })
  pct_of_genes_all$shortest_dist = shortest_dists
  pct_of_genes_all$pct_of_weights = do.call(rbind,pct_of_weights_list)[,1]

  #pct_of_genes_all$power = pct_of_genes_all$pct_of_genes * pct_of_genes_all$shortest_dist

  # make return object
  return_list =list(
    scores_per_leaf_region_all = scores_per_leaf_region_all,
    scores_per_target_level_region_all = scores_per_target_level_region_all,
    scores_per_voxel_annotated_list = scores_per_voxel_annotated_list,
    genes_per_voxel_list = genes_per_voxel_list,
    gene_set_power = pct_of_genes_all
  )
  return(return_list)

}

#' Annotate sets of cluster Markers with regions
#'
#' @param findRegion_result result list from findRegions_genesets
#' @param min_score min_score to filter for in findRegion_result enrichment (standard function return 0-1, and the default is set to 0.8)
#' @return  vector with all_children from nodes
#'
#' @export
#'
#'

summariseRegions_genesets = function(findRegion_result,min_score=0.8){

  all_gene_sets = rownames(findRegion_result$gene_set_power)
  result_region_list = list()
  for(gene_set in all_gene_sets){
    result_vec = vector()
    result_vec["gene_set"] = gene_set
    if(gene_set %in% colnames(findRegion_result$scores_per_leaf_region_all)){
      # if(max(findRegion_result$scores_per_target_level_region_all[,gene_set])[1]>min_score){
      # - Region (Target level)
      result_vec["Region"] = findRegion_result$scores_per_target_level_region_all[which(findRegion_result$scores_per_target_level_region_all[,gene_set] == max(findRegion_result$scores_per_target_level_region_all[,gene_set]))[1],2]
      # - Most likely region (all levels)
      result_vec["Region_specific"] = findRegion_result$scores_per_leaf_region_all[which(findRegion_result$scores_per_leaf_region_all[,gene_set] == max(findRegion_result$scores_per_leaf_region_all[,gene_set]))[1],2]
      result_vec["Region_specific_enrichment"] = findRegion_result$scores_per_leaf_region_all[which(findRegion_result$scores_per_leaf_region_all[,gene_set] == max(findRegion_result$scores_per_leaf_region_all[,gene_set]))[1],gene_set]
      # - Other likely regions (all levels, cutoff ?)
      top_5_scores = sort(findRegion_result$scores_per_leaf_region_all[,gene_set],decreasing = TRUE)[1:5]
      top_5_scores = top_5_scores[top_5_scores>min_score]
      top_5_scores = findRegion_result$scores_per_leaf_region_all[findRegion_result$scores_per_leaf_region_all[,gene_set] %in% top_5_scores,c(1,2,which(colnames(findRegion_result$scores_per_leaf_region_all)==gene_set))]
      top_5_scores = top_5_scores[order(top_5_scores[,gene_set],decreasing=TRUE),]
      result_vec["Region_specific_other"] = possible_regions = paste0(top_5_scores[-1,2],collapse = " | ")
      # }else{
      #   result_vec =c(result_vec,rep(NA,4))
      #   names(result_vec)[2:5] = c("Region","Region_specific","Region_specific_enrichment","Region_specific_other")
      # }
    }else{
      result_vec =c(result_vec,rep(NA,4))
      names(result_vec)[2:5] = c("Region","Region_specific","Region_specific_enrichment","Region_specific_other")
    }
    # - Power of markers used.
    result_vec["Region_confidence"] = findRegion_result$gene_set_power[gene_set,"pct_of_genes"]
    # result
    result_region_list[[gene_set]] = result_vec
  }
  result_region = as.data.frame(do.call(rbind,result_region_list))
  return(result_region)
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

#' Helper function for list intersection
#'
#' Adapted from : UpSetR::fromList but extended with rownames
#'
#' @param input list of gene sets (or any sets)
#' @return binary matrix of elements x sets
#'
#' @export
#'

intersect_from_list = function(input){
  # this code is adapted from : UpSetR::fromList
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  rownames(data) <- elements # added this to include rownames!
  return(data)
}



