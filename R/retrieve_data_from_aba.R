# # this script can be used to make a local copy of all allen brain ISH data
# # see also: https://github.com/AllenInstitute/cocoframer
#
# # target_folder
# target_folder = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/allen_brain/ish_data_voxel/"
# system(paste0("mkdir -p ",target_folder))
#
# # query_value
# query_value="energy" # what should be queried in cocoframer::get_aba_ish_data.  "energy", "density", "intensity", and "injection". See also: http://help.brain-map.org/display/mousebrain/API
#
# # get ids of available genes:
# exp_gene_relationships_df = cocoframer::get_exp_gene_relationships(plane = "coronal")
# ids_to_query = exp_gene_relationships_df$id[exp_gene_relationships_df$gene_symbol != ""]
# # save
# data.table::fwrite(exp_gene_relationships_df,file = paste0(target_folder,"aba_gene_to_id.tsv"),sep="\t",col.names = TRUE)
#
# # query all genes
# aba_ish_list = list()
# sttime=Sys.time()
# for(i in 1:length(ids_to_query)){
#   if(i %% 50 == 0){
#     message("Completed ",i-1," out of ",length(ids_to_query)," queries. Waiting for 2 seconds.")
#     Sys.sleep(time = 2)
#   }
#   current_id = ids_to_query[i]
#   temp <- tryCatch({
#     cocoframer::get_aba_ish_data(api_id=current_id, values = "energy")
#   },
#   error=function(cond) {
#     message("Query failed for:",current_id, " with error: ",cond)
#     return(NULL)
#   })
#   aba_ish_list[[current_id]] = temp
# }
# message("time taken: ",Sys.time()-sttime," ",attr(Sys.time()-sttime,"units"))
#
# # save list of 3d arrays as an R object
# saveRDS(aba_ish_list,paste0(target_folder,"aba_ish_list_",query_value,".rds"))
#
# ## make a matrix
# aba_ish_matrix = as.data.frame(do.call(cbind,aba_ish_list))
# data.table::fwrite(aba_ish_matrix,file = paste0(target_folder,"aba_ish_matrix_",query_value,".tsv"),sep="\t",row.names = TRUE, col.names = TRUE)
#
# ## retrieve other information
#
# ## need structures
# mba_ontology = cocoframer::get_mba_ontology()
# mba_ontology_flatten = cocoframer::flatten_mba_ontology(mba_ontology)
# mba_ontology_flatten$parent_structure_id = as.character(mba_ontology_flatten$parent_structure_id)
# mba_ontology_flatten$atlas_id = as.character(mba_ontology_flatten$atlas_id)
# mba_ontology_flatten$id = as.character(mba_ontology_flatten$id)
# data.table::fwrite(mba_ontology_flatten,file = paste0(target_folder,"mba_ontology_flatten.tsv"),sep="\t", col.names = TRUE)
#
# # get voxel structure information
# ccf_grid_annotation = cocoframer::get_ccf_grid_annotation()
# # save 3d array as an R object
# saveRDS(ccf_grid_annotation,paste0(target_folder,"aba_ccf_grid_annotation.rds"))
#
#
# ### make a matrix into coordinate system df for visualization and clustering
# # require(cubelyr)
# temp_3d_array = ccf_grid_annotation
# dimnames(temp_3d_array) <- list("x" = 1:dim(temp_3d_array)[1],
#                                 "y" = 1:dim(temp_3d_array)[2],
#                                 "z" = 1:dim(temp_3d_array)[3])
# temp_3d_cube <- cubelyr::as.tbl_cube(temp_3d_array)
# temp_2d_df <- dplyr::as_tibble(temp_3d_cube)
# temp_2d_df$structure_id = as.character(temp_2d_df$temp_3d_array)
# temp_2d_df$voxel_id = as.character(1:nrow(temp_2d_df))
# aba_matrix_annotation = dplyr::left_join(temp_2d_df[,c("voxel_id","x","y","z","structure_id")],mba_ontology_flatten[,c("id","name")],by=c("structure_id"="id"))
# # save
# data.table::fwrite(aba_matrix_annotation,file = paste0(target_folder,"aba_matrix_annotation.tsv"),sep="\t", col.names = TRUE)
#
#
#
