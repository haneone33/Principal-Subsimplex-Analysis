get_loading <- function(dendrogram.input){
  max_r = nrow(dendrogram.input$merges)
  modes = vector(mode = 'list', length = max_r)
  names(modes) = rownames(dendrogram.input$merges)
  
  for(r in (max_r-1):1){
    v1 = dendrogram.input$vertices[[paste0('r=',r+1)]][dendrogram.input$merges[paste0('r=',r),'v1'],]
    v2 = dendrogram.input$vertices[[paste0('r=',r+1)]][dendrogram.input$merges[paste0('r=',r),'v2'],]
    modes[[paste0('r=',r)]] = list(v1 = v1, v2 = v2)
  }
  
  return(modes)
}