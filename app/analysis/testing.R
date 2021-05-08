significant.velocities <- function(df_vt.smooth) {
    wilcox.test(df_vt.smooth$vel.ma ~ df_vt.smooth$group, conf.int = TRUE)
}

chisq.multicomp <- function(x,p.method="fdr") {
  x <- sort(x)
  fun.p <- function(i,j) {
    xi <- x[i]
    xj <- x[j]
    suppressWarnings(chisq.test(c(xi, xj)))$p.value
  }
  tab.p <- pairwise.table(fun.p,as.character(x),p.adjust.method=p.method)
  call <- match.call()
  dname.x <- if(length(call$x)==1) {call$x} else {paste(call$x[1],"(",paste(call$x[-1],collapse=","),")",sep="")}
  result <- list(method="chi-squared tests",data.name=dname.x,p.adjust.method=p.method,p.value=tab.p)
  class(result) <- "pairwise.htest"
  return(result)
}

segmentation.significance <- function(segments, method, p.correct) {
    ngroups <- segments %>% select(group) %>% distinct()
    if (nrow(ngroups) < 2) {
        return(data.frame(x = "Could not quantify significance. Control group needed."))
    }

    df.seg <- segments %>% group_by(group) %>% count(cluster) %>% mutate(total = sum(n), prop = n / sum(n))

    mod.1 <-
        glm(group ~ cluster * begin * end,
            data = segments,
            family = "binomial")
    mod.2 <-
        glm(group ~ begin * end, data = segments, family = "binomial")
    mod.3 <-
        glm(group ~ cluster, data = segments, family = "binomial")
    mod.5 <-
        glm(group ~ cluster * prop, data = df.seg, family = "binomial")
    mod.6 <-
        glm(group ~ cluster, data = df.seg, family = "binomial")
    mod.4 <- glm(group ~ 1, data = segments, family = "binomial")

    print(anova(mod.1, mod.2, test = method))
    print(anova(mod.1, mod.4, test = method))
    print(anova(mod.2, mod.4, test = method))
    print(anova(mod.5, mod.6, test = method))
    anova(mod.3, mod.4, test = method)
}

spectrum.significance <-
    function(df_freqs,
             range) {
        df_freqs <- df_freqs %>% mutate(spec.f = spec.f / max(spec.f)) %>%
            filter(spec.f < range[2], spec.f > range[1])

        m1.r <-
            glm(factor(group) ~ spec.f * spec.s ,
                data = df_freqs,
                family = "binomial")
        m2.r <-
            glm(factor(group) ~ spec.f + spec.s ,
                data = df_freqs,
                family = "binomial")
        m3.r <-
            glm(factor(group) ~ 1 ,
                data = df_freqs,
                family = "binomial")
        anova(m1.r, m2.r, m3.r, test = "Chisq")
    }