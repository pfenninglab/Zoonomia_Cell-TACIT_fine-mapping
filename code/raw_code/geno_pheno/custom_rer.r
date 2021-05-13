#These are functions from the RERconverge package that are very lightly modified to suit my use case.

library(RERconverge)

getAllCladeEdges=function(tree, AncEdge){
  node=tree$edge[AncEdge,2]
  #get descendants
  iid=getDescendants(tree, node)
  #find their edges
  iim=match(iid, tree$edge[,2])
  iim
}

#This function is unchanged, but duplicated because it is tagged as an internal function.
#' @keywords internal
nameEdges=function(tree){
  nn=character(nrow(tree$edge))
  iim=match(1:length(tree$tip.label), tree$edge[,2])
  nn[iim]=tree$tip.label
  nn
}

#This function is unchanged, but duplicated because it is tagged as an internal function.
#' @keywords internal
inferUnidirectionalForegroundClades <- function(tree, fgd = NULL, ancestralOnly = F){
  finaltree <- tree
  finaltree$edge.length <- rep(0, length(tree$edge.length))
  finaltree$edge.length[nameEdges(finaltree) %in% fgd] <- 1
  #figure out node depth - terminal nodes have depth of 1; higher numbers indicate ancestral nodes;
  nodedepths <- node.depth(finaltree)
  edgeterminalnodedepths <- nodedepths[finaltree$edge[,2]]
  #going from 1-away from terminal ancestral branch to the base of the tree, figure out branches where all downstream clades are foreground
  for(inode in sort(unique(edgeterminalnodedepths))[-1]){
    edgesToDo <- which(edgeterminalnodedepths == inode)
    for(edgeindex in edgesToDo){
      clade.edges = getAllCladeEdges(finaltree, edgeindex)
      if(all(finaltree$edge.length[clade.edges]==1)){
        finaltree$edge.length[edgeindex] <- 1
      }
    }
  }
  if(ancestralOnly){
    for(edgeii in 1:length(finaltree$edge.length)){
      if(finaltree$edge.length[edgeii] == 1){
        if(nameEdges(finaltree)[edgeii]==""){
          clade.edges = setdiff(getAllCladeEdges(finaltree, edgeii), edgeii)
          finaltree$edge.length[clade.edges] <- 0
        }
      }
    }
  }
  finaltree
}

#This function is changed to take a single phylo tree as opposed to a treesObj
#' Creates a binary trait tree from a set of foreground species.
#' @param foreground. A character vector containing the foreground species
#' @param tree A tree
#' @param collapse2anc Put all the weight on the ancestral branch when the trait appears on a whole clade
#' (redundant to "clade", kept for backwards compatibility)
#' @param plotTree Plot a tree representation of the result
#' @param wholeClade Whether to implement the weighted edge option across
#' all members of a foreground clade (redundant to "clade", kept for backwards compatibility)
#' @param clade A character string indicating which branches within the clade
#' containing the foreground species should be set to foreground. Must be one
#' of the strings "ancestral", "terminal", "all".
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @param weighted if set to TRUE weights foreground edges belonging to the same clade such that their branch lengths sum up to 1 (only done for clade options "all" and "terminal").
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A tree with edge.lengths representing phenotypic states
#' @export
customForeground2Tree = function(foreground,tree, plotTree=F, clade=c("ancestral","terminal","all"), weighted = F, transition = "unidirectional", useSpecies=NULL){
  clade <- match.arg(clade)
  res = tree #treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
                                                                                   collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  } else {
    useSpecies = res$tip.label
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0,length(res$edge.length))
  if(clade == "terminal"){
    res$edge.length[nameEdges(res) %in% foreground] = 1
    names(res$edge.length) = nameEdges(res)
  }else if(clade == 'ancestral'){
    weighted = F
    if(transition == 'bidirectional'){
      res <- inferBidirectionalForegroundClades(res,foreground,ancestralOnly =T)
    }else{
      res <- inferUnidirectionalForegroundClades(res,foreground,ancestralOnly =T)
    }
  }else{
    if(transition == 'bidirectional'){
      res <- inferBidirectionalForegroundClades(res,foreground,ancestralOnly = F)
    }else{
      res <- inferUnidirectionalForegroundClades(res,foreground,ancestralOnly = F)
    }
  }
  if(weighted){
    if(clade == 'all'){
      tobeweighted <- rep(TRUE,length(res$edge.length))
      tobeweighted[res$edge.length == 0] <- FALSE
      while(sum(tobeweighted)>0){
        edgetodo <- which(tobeweighted == T)[1]
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if(length(clade.down.edges) > 1){
          clade.edges = c(clade.down.edges, edgetodo)
          clade.edges.toweight <- clade.edges[res$edge.length[clade.edges] == 1]
          res$edge.length[clade.edges.toweight] <- 1.0/(length(clade.edges.toweight))
          tobeweighted[clade.edges] <- FALSE
        } else{
          tobeweighted[clade.down.edges] <- FALSE
        }
      }
    } else if(clade == 'terminal'){
      tobeweightededgeterminalnode <- unique(res$edge[(res$edge[,2] %in% c(1:length(res$tip.label))),1])
      tobeweighted <- setdiff(match(tobeweightededgeterminalnode, res$edge[,2]), NA)
      for(edgetodo in tobeweighted){
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if(all(res$edge.length[clade.down.edges]==1)){
          res$edge.length[clade.down.edges] <- 0.5
        }
      }
    }
  }
  if(plotTree){
    res2=res
    mm=min(res2$edge.length[res2$edge.length>0])
    res2$edge.length[res2$edge.length==0]=max(0.02,mm/20)
    plot(res2, main = paste0("Clade: ",clade,'\nTransition: ',transition,'\nWeighted: ',weighted), cex = 0.5)
    if(weighted){
      labs <- round(res$edge.length,3)
      labs[labs == 0] <- NA
      edgelabels(labs, col = 'black', bg = 'transparent', adj = c(0.5,-0.5),cex = 0.4,frame='n')
    }
  }
  res
}

#This function is changed to take a single phylo tree as opposed to a treesObj, and to change how counts of foreground tip and internal nodes are passed.
#'Generates a permulated phenotype vector whose phylogeny matches a desired structure.  User may specify the number of foreground branches that are internal branches.
#' @param tree A tree
#' @param root Species on which to root the master tree
#' @param phenvec Named vector of 1's and 0's representing phenotype values for each species
#' @param fgnum Total number of foreground species (including ancestors) - only required if internal foreground branches are required
#' @return A vector of permulated foreground species
#' @export
customSimBinPhenoVec=function(tree, root, phenvec, fgnum=NULL){
  blsum=0
  tips = sum(phenvec)
  if(is.null(fgnum)){
    fgnum=tips
  }
  internal=fgnum-tips
  while(blsum!=fgnum){
    t=tree #root.phylo(tree, root, resolve.root = T) #For now, let's use the ancestral root
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    t=customForeground2Tree(top, tree, clade="all", plotTree = F)
    blsum=sum(t$edge.length)
  }
  # plot(t)
  return(top)
}
