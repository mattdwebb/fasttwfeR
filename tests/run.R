args <- commandArgs(TRUE)
root <- if (length(args)) args[1L] else "."
if (file.exists(file.path(root, "R", "estimator.R"))) {
  source(file.path(root, "R", "estimator.R"))
} else {
  ns <- asNamespace("fasttwfeR")
  for (name in c(".run", ".prepare", ".assignments", ".family"))
    assign(name, get(name, envir = ns))
}
fixtures <- file.path(root, "tests", "fixtures")
if (!dir.exists(fixtures)) fixtures <- "fixtures"
f <- .run
d <- read.csv(file.path(fixtures, "data.csv"))
aa <- as.matrix(read.csv(file.path(fixtures, "assignments.csv"), header = FALSE))
checks <- 0L
check <- function(value) { stopifnot(value); checks <<- checks + 1L }
maxdiff <- 0
near <- function(x,y,tol=1e-9) {
  delta <- max(abs(x-y))
  maxdiff <<- max(maxdiff, delta)
  check(delta < tol)
}
fails <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error=function(e) conditionMessage(e))
  check(!is.null(e) && grepl(pattern,e))
}
for (cov in list(character(), "x")) {
  suffix <- if(length(cov)) "_x" else ""
  expected <- read.csv(file.path(fixtures,paste0("expected",suffix,".csv")))
  dels <- read.csv(file.path(fixtures,paste0("deletions",suffix,".csv")))$deleted
  rr <- as.matrix(read.csv(file.path(fixtures,paste0("randomized",suffix,".csv"))))
  fit <- f(d,"y","time","state",gvar="g",covariates=cov,
           inference=c("native","jackknife","ri_beta","ri_t_native"),assignments=aa)
  for (name in names(expected)) near(fit[[name]],expected[[name]])
  near(fit$deletions,dels)
  near(fit$randomized_statistics,rr)
  observed <- c(fit$estimate,fit$native_stat)
  for(j in seq_len(2)) near(fit[[paste0(colnames(rr)[j],"_p")]],
    (1+sum(abs(rr[,j])>=abs(observed[j])-1e-12))/(nrow(rr)+1))
  fit2 <- f(d[nrow(d):1,],"y","time","state",treatment="d",covariates=cov,inference="point")
  near(fit$estimate,fit2$estimate)
  check(fit$jackknife_df==7)
  y2 <- d; y2$y <- 2*d$y+sin(d$x)
  a1 <- f(y2,"y","time","state",gvar="g",covariates=cov,inference="point")$estimate
  y2$y <- sin(d$x)
  a2 <- f(y2,"y","time","state",gvar="g",covariates=cov,inference="point")$estimate
  near(a1,2*fit$estimate+a2)
}
bad <- d; bad$d[1] <- 1
fails(f(bad,"y","time","state",treatment="d"),"inconsistent")
bad <- d; bad$y[1] <- NA
fails(f(bad,"y","time","state",gvar="g"),"missing")
fails(f(d,"y","time","state",gvar="g",inference="madeup"),"UNSUPPORTED")
fails(f(d,"y","time","state",gvar="g",inference="ri_beta",assignments=rbind(aa,aa[1,])),"duplicate")
fails(f(d,"y","time","state",gvar="g",inference="ri_beta",assignments=matrix(c(0,0,0,3,3,4,4,0),1)),"observed")
fails(f(d,"y","time","state",gvar="g",inference="ri_beta",assignments=aa,null=1),"UNSUPPORTED")
fails(f(d,"y","time","state",gvar="g",data_mode="panel"),"requires unit")
bad <- d; bad$x2 <- bad$x
if (.family != "flex") fails(f(bad,"y","time","state",gvar="g",covariates=c("x","x2")),"rank")
a <- .prepare(d,"y","time","state",NULL,"g",character(),NULL,NULL,"rcs")$a
set.seed(821); before <- .Random.seed
draws <- .assignments(a,3,92)
check(identical(before,.Random.seed))
check(identical(draws,.assignments(a,3,92)))
check(nrow(.assignments(c(0,0,3,3),999,1))==5)
bad <- d[d$state!=4,]
surviving <- f(bad,"y","time","state",gvar="g",inference="jackknife")
physical <- vapply(sort(unique(bad$state)),function(i)
  f(bad[bad$state!=i,],"y","time","state",gvar="g",inference="point")$estimate,numeric(1))
near(surviving$deletions,physical)
if(.family %in% c("twfe","did2s")) {
  wfit <- f(d,"y","time","state",gvar="g",covariates="x",weights="w",inference="native")
  X <- model.matrix(~factor(state)+factor(time)+x,data=d)
  if(.family=="twfe") {
    m <- lm(y~factor(state)+factor(time)+x+d,data=d,weights=w)
    near(wfit$estimate,unname(coef(m)["d"]))
    Z <- model.matrix(m)
    direction <- solve(crossprod(Z,d$w*Z),diag(ncol(Z))[,ncol(Z)])
    scores <- rowsum(d$w * residuals(m)*as.vector(Z%*%direction),d$state)
    se <- sqrt(sum(scores^2)*8/7*(nrow(d)-1)/(nrow(d)-ncol(Z)))
    near(wfit$native_se,se)
    panel <- f(d,"y","time","state",gvar="g",covariates="x",unit="unit",inference="native")
    m <- lm(y~factor(unit)+factor(time)+x+d,data=d)
    near(panel$estimate,unname(coef(m)["d"]))
  } else {
    m <- lm(y~factor(state)+factor(time)+x,data=d[d$d==0,],weights=w)
    near(wfit$estimate,weighted.mean(d$y[d$d==1]-predict(m,d[d$d==1,]),d$w[d$d==1]))
  }
} else fails(f(d,"y","time","state",gvar="g",weights="w",inference="native"),"UNSUPPORTED")
cat(.family,":",checks,"checks passed; max discrepancy",maxdiff,"\n")
