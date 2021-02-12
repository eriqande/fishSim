#' partially-parallelized findRelatives() hacked by ECA to only go out to the great-great grandparents
#'
#' findRelativesPar is a partially-parallelized version of findRelatives. Its use is
#' exactly the same, but it requires libraries 'foreach', 'parallel', and 'doParallel'.
#' 'findRelativesAlt' uses the sample indicator in indiv[,9] to decide which individuals
#' to compare, whereas 'findRelativesPar' compares between all individuals in 'sampled'.
#'
#' Lookup operations to find each member's ancestors are parallelized, but comparisons
#' between ancestor-sets are not. On a test-set of 100 sampled individuals, this partial-
#' parallelization reduced runtime from 55 seconds to 22 seconds.
#' findRelatives takes a set of 'sampled' individuals and a population simulation output
#' (i.e., an 'indiv' object, of the kind output by mort() ). It returns a data.frame
#' for each pair of sampled individuals, with columns:
#' 1) Var1: the first individual's ID
#' 2) Var2: the second individual's ID
#' 3) related: TRUE for each pair that with one or more shared ancestors within seven
#'    ancestral generations (i.e., great-great-great-great grandparents), otherwise FALSE
#' 4) totalRelatives: numeric indicator of the total number of shared ancestors found
#' 5) OneTwo, OneThree, ThreeFour, etc.: Numeric values, indicating the number of shared
#'    ancestors by relationship class. For instance, a shared relative in the OneTwo class
#'    indicates that one member of the pair is the other member's parent, a shared relative
#'    in the OneThree class indicates that one member is the other's grandparent, and a shared
#'    in the ThreeFour class indicates that one individual's grandparent is the other's great-
#'    grandparent. If a pair shares an ancestor in the TwoThree class, they necessarily also
#'    share two ancestors in the ThreeFour class and four ancestors in the FourFive class,
#'    and so on. Note that relationship classes are identical by reversal - ThreeFour is the
#'    same as FourThree, and so only relationship classes with increasing order are presented
#'    (i.e., ThreeFour and OneFive are in the output, but not ThreeTwo).
#' Note that, if there is an object called 'ancestors' in the global environment, the
#' foreach() loops may refer to that copy of 'ancestors', rather than the one generated inside
#' findRelativesPar(). This is a known bug. Rename and delete 'ancestors'.
#' If it occurs, this bug will cause the following error:
#' Error in cbind(ancestors, parents.o, grandparents.o, ggrandparents.o:
#' number of rows of matrices must match (see arg 3).
#' @param indiv A matrix of individuals, as from mort(), but which will need to contain
#'              long strings of parent-offspring relationships, so will most likely come
#'              from a multi-generation simulation.
#' @param sampled TRUE or FALSE. If TRUE, compares only individuals marked as sampled in
#'                indiv[,9]. If FALSE, compares all individuals in 'indiv'.
#' @param verbose TRUE or FALSE. If TRUE, prints a table of sampling
#'                years for sampled individuals.
#' @param nCores the number of cores to use for parallel processes. Defaults to one less than
#'               the number of cores on the machine.
#' @param delimitIndiv TRUE/FALSE. Lookups can be sped up markedly by first delimiting
#'                     indiv to only those animals that exist as parents, or that are
#'                     marked as sampled. Default TRUE.
#' @seealso [fishSim::findRelatives()]
#' @seealso [fishSim::capture()]
#' @export
findRel4gen <- function(indiv, sampled = TRUE, verbose = TRUE, nCores = detectCores()-1,
                             delimitIndiv = TRUE) {

  registerDoParallel(nCores)

  if (sampled) {
    if (sum(!is.na(indiv[,9])) == 0)
      stop("no sampled individuals")
    if (verbose)
      print(data.frame(table(indiv[!is.na(indiv[,9]),9], dnn = "Sample Year")))
    sampled <- indiv[!is.na(indiv[,9]),1]

  } else sampled <- indiv[,1]

  if(delimitIndiv) {
    keepers <- indiv$Me %in% sampled | indiv$Me %in% indiv$Mum | indiv$Me %in% indiv$Dad
    indiv <- indiv[keepers,]
  }

  ancestors <- matrix(data = sampled, nrow = length(sampled))

  parents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% {
    parents(ancestors[i,1], indiv)
  }
  print(paste("parents found at ", Sys.time(), sep = ""))
  grandparents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% {
    grandparents(ancestors[i,1], indiv)
  }
  print(paste("grandparents found at ", Sys.time()))
  ggrandparents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% {
    great.grandparents(ancestors[i,1], indiv)
  }
  print(paste("great-grandparents found at ", Sys.time(), sep = ""))
  gggrandparents.o <- foreach(i = 1:nrow(ancestors), .combine = rbind) %dopar% {
    great2.grandparents(ancestors[i,1], indiv)
  }

  print(paste("great-great-grandparents found at ", Sys.time(), sep = ""))


  ancestors <- cbind(ancestors, parents.o, grandparents.o, ggrandparents.o, gggrandparents.o)
  colnames(ancestors) <- c("self","father", "mother", # self and parents
                           "FF", "FM", "MF", "MM", # grandparents
                           "FFF","FFM","FMF","FMM","MFF","MFM","MMF","MMM", #great-grandparents
                           "FFFF","FFFM","FFMF","FFMM","FMFF","FMFM","FMMF","FMMM",
                           "MFFF","MFFM","MFMF","MFMM","MMFF","MMFM","MMMF","MMMM" #gg-grandparents
                           )

  expand.grid.unique <- function(x, y, include.equals=FALSE) {
    x <- unique(x)
    y <- unique(y)
    g <- function(i) {
      z <- setdiff(y, x[seq_len(i-include.equals)])
      if(length(z)) cbind(x[i], z, deparse.level=0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  } ## with thanks to stack overflow user Ferdinand.kraft

  pairs <- expand.grid.unique(ancestors[,1], ancestors[,1])
  colnames(pairs) <- c("Var1", "Var2")

  ##    pairs <- expand.grid(ancestors[,1], ancestors[,1])
  ##    pairs <- pairs[pairs$Var1 != pairs$Var2,] ## remove self-comparisons

  related <- c(rep(NA, nrow(pairs)))
  totalRelatives <- c(rep(NA, nrow(pairs))) ##mainly here for imagined post-hoc diagnostics

  OneTwo <- c(rep(NA, nrow(pairs)))
  OneThree <- c(rep(NA, nrow(pairs)))
  OneFour <- c(rep(NA, nrow(pairs)))
  OneFive <- c(rep(NA, nrow(pairs)))
  TwoTwo <- c(rep(NA, nrow(pairs)))
  TwoThree <- c(rep(NA, nrow(pairs)))
  TwoFour <- c(rep(NA, nrow(pairs)))
  TwoFive <- c(rep(NA, nrow(pairs)))
  ThreeThree <- c(rep(NA, nrow(pairs)))
  ThreeFour <- c(rep(NA, nrow(pairs)))
  ThreeFive <- c(rep(NA, nrow(pairs)))
  FourFour <- c(rep(NA, nrow(pairs)))
  FourFive <- c(rep(NA, nrow(pairs))) ## PC
  FiveFive <- c(rep(NA, nrow(pairs)))  ## PC

  for(i in 1:length(related)) {
    allAncestors <- ancestors[ancestors[,1] == pairs[i,1],1:31] %in%
      ancestors[ancestors[,1] == pairs[i,2],1:31]
    indi1 <- pairs[i,1] ## 1st individual's self
    indi2 <- pairs[i,2] ## 2nd individual's self
    indi1Par <- ancestors[ancestors[,1] == pairs[i,1],2:3] ## 1st indiv's parents
    indi2Par <- ancestors[ancestors[,1] == pairs[i,2],2:3] ## 2nd indiv's parents
    indi1GP <- ancestors[ancestors[,1] == pairs[i,1],4:7] ## 1st indiv's grandparents
    indi2GP <- ancestors[ancestors[,1] == pairs[i,2],4:7] ## 2nd indiv's grandparents
    indi1GGP <- ancestors[ancestors[,1] == pairs[i,1],8:15] ## 1st indiv's g-grandparents
    indi2GGP <- ancestors[ancestors[,1] == pairs[i,2],8:15] ## 2nd indiv's g-grandparents
    indi1GGGP <- ancestors[ancestors[,1] == pairs[i,1],16:31] ## 1st indiv's gg-grandparents
    indi2GGGP <- ancestors[ancestors[,1] == pairs[i,2],16:31] ## 2nd indiv's gg-grandparents
#    indi1GGGGP <- ancestors[ancestors[,1] == pairs[i,1],32:63] ## 1st indiv's ggg-grandparents
#    indi2GGGGP <- ancestors[ancestors[,1] == pairs[i,2],32:63] ## 2nd indiv's ggg-grandparents
#    indi1GGGGGP <- ancestors[ancestors[,1] == pairs[i,1],64:127] ## 1st indiv's gggg-grandparents
#    indi2GGGGGP <- ancestors[ancestors[,1] == pairs[i,2],64:127] ## 2nd indiv's gggg-grandparents

    related[i] <- any(allAncestors)
    totalRelatives[i] <- sum(allAncestors)

    ## Identify pairs where one indiv is the other indiv's n-parent
    OneTwo[i] <- sum(c(indi1 %in% indi2Par, indi2 %in% indi1Par))
    OneThree[i] <- sum(c(indi1 %in% indi2GP, indi2 %in% indi1GP))
    OneFour[i] <- sum(c(indi1 %in% indi2GGP, indi2 %in% indi1GGP))
    OneFive[i] <- sum(c(indi1 %in% indi2GGGP, indi2 %in% indi1GGGP))
#    OneSix[i] <- sum(c(indi1 %in% indi2GGGGP, indi2 %in% indi1GGGGP))
#    OneSeven[i] <- sum(c(indi1 %in% indi2GGGGGP, indi2 %in% indi1GGGGGP))

    ## Identify pairs where one indiv's parent is the other's n-parent
    TwoTwo[i] <- sum(c(indi1Par %in% indi2Par))
    TwoThree[i] <- sum(c(indi1Par %in% indi2GP, indi2Par %in% indi1GP))
    TwoFour[i] <- sum(c(indi1Par %in% indi2GGP, indi2Par %in% indi1GGP))
    TwoFive[i] <- sum(c(indi1Par %in% indi2GGGP, indi2Par %in% indi1GGGP))
#    TwoSix[i] <- sum(c(indi1Par %in% indi2GGGGP, indi2Par %in% indi1GGGGP))
#    TwoSeven[i] <- sum(c(indi1Par %in% indi2GGGGGP, indi2Par %in% indi1GGGGGP)) ## PC

    ## Identify pairs where one indiv's grandparent is the other's n-parent
    ThreeThree[i] <- sum(c(indi1GP %in% indi2GP))
    ThreeFour[i] <- sum(c(indi1GP %in% indi2GGP, indi2GP %in% indi1GGP))
    ThreeFive[i] <- sum(c(indi1GP %in% indi2GGGP, indi2GP %in% indi1GGGP))
#    ThreeSix[i] <- sum(c(indi1GP %in% indi2GGGGP, indi2GP %in% indi1GGGGP)) ## PC
#    ThreeSeven[i] <- sum(c(indi1GP %in% indi2GGGGGP, indi2GP %in% indi1GGGGGP)) ## PC

    ## Identify pairs where one indiv's g-grandparent is the other's n-parent
    FourFour[i] <- sum(c(indi1GGP %in% indi2GGP))
    FourFive[i] <- sum(c(indi1GGP %in% indi2GGGP, indi2GGP %in% indi1GGGP)) ## PC
#    FourSix[i] <- sum(c(indi1GGP %in% indi2GGGGP, indi2GGP %in% indi1GGGGP)) ## PC
#    FourSeven[i] <- sum(c(indi1GGP %in% indi2GGGGGP, indi2GGP %in% indi1GGGGGP)) ## PC

    ## Identify pairs where one indiv's gg-grandparent is the other's n-parent
    FiveFive[i] <- sum(c(indi1GGGP %in% indi2GGGP)) ## PC
#    FiveSix[i] <- sum(c(indi1GGGP %in% indi2GGGGP, indi2GGGP %in% indi1GGGGP)) ## PC
#    FiveSeven[i] <- sum(c(indi1GGGP %in% indi2GGGGGP, indi2GGGP %in% indi1GGGGGP)) ## PC

    ## Identify pairs where one indiv's ggg-grandparent is the other's n-parent
#    SixSix[i] <- sum(c(indi1GGGGP %in% indi2GGGGP))  ## PC
#    SixSeven[i] <- sum(c(indi1GGGGP %in% indi2GGGGGP, indi2GGGGP %in% indi1GGGGGP)) ## PC

    ## Identify pairs where one indiv's gggg-grandparent is the other's gggg-grandparent
#    SevenSeven[i] <- sum(c(indi1GGGGGP %in% indi2GGGGGP)) ## possibly cull
    if(i %% 1000 == 0) {
      cat( '\r', i, " of ", length(related), " comparisons", sep = ""); flush.console()
    }
  }

  pairs <- data.frame(pairs, related, totalRelatives,
                      OneTwo, OneThree, OneFour, OneFive, #OneSix, OneSeven,
                      TwoTwo, TwoThree, TwoFour, TwoFive, #TwoSix,
                      #TwoSeven,  ## possibly cull
                      ThreeThree, ThreeFour, ThreeFive,
                      #ThreeSix, ThreeSeven, ## possibly cull
                      FourFour,
                      FourFive, #FourSix, FourSeven, ## possibly cull
                      FiveFive #FiveSix, FiveSeven, ## possibly cull
                      #SixSix, SixSeven, SevenSeven  ## possibly cull
  ) ## there are a lot of identified relative-classes
  ## with absolutely minimal shared DNA here, and they are
  ## also the slowest to compute. 'possibly cull' (or 'PC') indicates things that could be
  ## removed for a code speed-up (measured at about 10% of runtime, on a
  ## 100-individual sample), at potentially not much inferential cost.

  stopImplicitCluster()  ## cleanup the cluster afterwards; not important but clean.
  return(pairs)
}

