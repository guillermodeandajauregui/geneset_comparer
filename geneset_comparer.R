
deregulation_comparer<-function(x){
#takes a list of sets as character vectors (for instance, deregulated genes in phenotypes)
#returns a list of lists of sets as character vectors representing their position on a Venn diagram for as many sets
#designed to compare 4 sets. 
  #function for intersections beyond "three out of four" are not implemented. 
  
all_deregulated <- Reduce(union, x)
only_in_subtype <-unintersected(x)
in_two_subtypes<-unique_intersector(x)
three_out_of_four<-intersect_all_but_one(x)
core_deregulated_genes <- Reduce(intersect, x)

resultados<-list(all_deregulated = all_deregulated, only_in_subtype = only_in_subtype, in_two_subtypes = in_two_subtypes, three_out_of_four = three_out_of_four, core_deregulated_genes = core_deregulated_genes)

return(resultados)
}

#functions needed:
#unintersected 
#unique_intersectors
#intersect_all_but_one
##
unintersected<-function(x){
  
  #x is a list
  setunion<-Reduce(union, x)
  listy <- list()
  
  #elements that do not intersect for each set 
  for (i in seq_along(x)) {
    ex_minus_i <- x[-i]
    union_ex_minus_i <- Reduce(union, ex_minus_i)
    listy[[i]] <- setdiff(setunion, union_ex_minus_i)
  }
  names(listy)<-paste0("ONLY_IN_", names(x))
  return(listy)
}
##
unique_intersector<-function(x){
  #function. takes a named list of character vectors (genesets). 
  #Returns a named list of intersections between all pairs of genesets as character vectors.
  lista_intersecciones<-list()
  k <- combn(x, 2)
  namek<-combn(names(x),2)
  
  for (i in seq_along(k[1,])){
    lista_intersectors <- x[namek[,i]]
    lista_outersectors <- x[names(x)%in%names(lista_intersectors) == FALSE]
    
    intersect_intersectors <- Reduce(intersect, lista_intersectors)
    union_outersectors <- Reduce(union, lista_outersectors)
    lista_intersecciones[[ paste0(namek[,i][1], "_", namek[,i][2], "_UNIQUE_INTERSECTION") ]] <- setdiff(intersect_intersectors, union_outersectors)
    
  }
  return(lista_intersecciones)
}
##
intersect_all_but_one<-function(x){
  
  #x is a list
  setunion<-Reduce(union, x)
  listy <- list()
  
  #elements that do not intersect for each set 
  for (i in seq_along(x)) {
    ex_minus_i <- x[-i]
    q <- Reduce(intersect, ex_minus_i)
    listy[[i]] <-setdiff(q, x[[i]])
  }
  names(listy)<-paste0("ALL_BUT_", names(x))
  return(listy)
}
