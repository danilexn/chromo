find.pc.causality <- function(X, pval, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    updateProgress("Step 1/3: Computing PC algorithm")
  }
  stuffStat <- list(C = cor(as.matrix(X)), n = nrow(as.matrix(X)))
  pc.fit <- pc(stuffStat, indepTest = gaussCItest, p = ncol(as.matrix(X)), alpha = pval)
  t(as(pc.fit, "matrix"))
}

find.vlgc.causality <- function(X, nlags, pval, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    updateProgress("Step 3/3: Computing VL-GC algorithm")
  }
  out <- multipleVLGrangerFunc(as.matrix(X), maxLag = as.numeric(nlags), alpha = pval)
  out$adjMat
}

find.vlte.causality <- function(X, nlags, pval, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    updateProgress("Step 2/3: Computing VL-TE algorithm")
  }
  out <- multipleVLTransferEntropy(as.matrix(X), maxLag = as.numeric(nlags), autoLagflag = TRUE, alpha = pval)
  out$adjMat
}

causality.global <- function(df, vars, pval, nlags, p.adj.method, sel.group, progress = NULL){
    # Mutate variables of signal and morphology
    df_filter <- df[,c("particle", "group", vars)] %>% filter(group == sel.group) %>% select(-group)
    nparts <- nrow(df_filter %>% dplyr::select(particle) %>% distinct())
    pval <- p.adjust.inv(pval, p.adj.method, nparts)
    updateProgress <- function(detail = NULL) {
      progress$inc(amount = 1/(nparts*3), detail = detail)
    }

    # Case 1: PC causality
    list.causality.pc <- df_filter %>% dplyr::group_by(particle) %>%
        group_map( ~ find.pc.causality(., pval, updateProgress))

    adj.mat.pc <- round(Reduce('+', list.causality.pc) / nparts, 2)

    # Case 2: VLTE causality
    list.causality.vlte <- df_filter %>% dplyr::group_by(particle) %>%
        group_map( ~ find.vlte.causality(., nlags, pval, updateProgress))

    adj.mat.vlte <- round(Reduce('+', list.causality.vlte) / nparts, 2)

    # Case 3: VLGC causality
    list.causality.vlgc <- df_filter %>% dplyr::group_by(particle) %>%
        group_map( ~ find.vlgc.causality(., nlags, pval, updateProgress))

    adj.mat.vlgc <- round(Reduce('+', list.causality.vlgc) / nparts, 2)

    return(list(adj.mat.pc, adj.mat.vlte, adj.mat.vlgc))
}