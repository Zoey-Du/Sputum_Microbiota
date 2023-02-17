library(ape)
library(lattice)
library(matrixStats)
library(permutes)
library(lmPerm)
library(permute)
library(vegan)
library(GUniFrac)
adonis2N <- function (formula, data, permutations = 999, method = "bray", tree = NULL,
          sqrt.dist = FALSE, add = FALSE, by = "terms", parallel = getOption("mc.cores"), 
          ...) 
{
  if (!is.null(by)) 
    by <- match.arg(by, c("terms", "margin"))
  YVAR <- formula[[2]]
  lhs <- eval(YVAR, environment(formula), globalenv())
  environment(formula) <- environment()
  if ((is.matrix(lhs) || is.data.frame(lhs)) && isSymmetric(unname(as.matrix(lhs)))) 
    lhs <- as.dist(lhs)
  if (!inherits(lhs, "dist")) 
    lhs <- vegdist2(as.matrix(lhs), method = method, tree = tree,...)
  if (sqrt.dist) 
    lhs <- sqrt(lhs)
  if (is.logical(add) && isTRUE(add)) 
    add <- "lingoes"
  if (is.character(add)) {
    add <- match.arg(add, c("lingoes", "cailliez"))
    if (add == "lingoes") {
      ac <- addLingoes(as.matrix(lhs))
      lhs <- sqrt(lhs^2 + 2 * ac)
    }
    else if (add == "cailliez") {
      ac <- addCailliez(as.matrix(lhs))
      lhs <- lhs + ac
    }
  }
  if (!missing(data)) 
    formula <- terms(formula, data = data)
  formula <- update(formula, lhs ~ .)
  if (missing(data)) 
    data <- model.frame(delete.response(terms(formula)))
  sol <- vegan:::adonis0(formula, data = data, method = method)
  out <- anova(sol, permutations = permutations, by = by, parallel = parallel)
  att <- attributes(out)
  out <- rbind(out, Total = c(nobs(sol) - 1, sol$tot.chi, NA, 
                              NA))
  out <- cbind(out[, 1:2], R2 = out[, 2]/sol$tot.chi, out[, 
                                                          3:4])
  att$heading[2] <- deparse(match.call(), width.cutoff = 500L)
  att$names <- names(out)
  att$row.names <- rownames(out)
  attributes(out) <- att
  out
}



vegdist2 <- function (x, method = "bray", tree = NULL, binary = FALSE, diag = FALSE, upper = FALSE, 
                      na.rm = FALSE, ...) 
{
  if (method != 17 & method != 18) x1 = x
  ZAP <- 1e-15
  if (!is.na(pmatch(method, "euclidian"))) 
    method <- "euclidean"
  METHODS <- c("manhattan", "euclidean", "canberra", "bray", 
               "kulczynski", "gower", "morisita", "horn", "mountford", 
               "jaccard", "raup", "binomial", "chao", "altGower", "cao", 
               "mahalanobis", "dw", "du")
  method <- pmatch(method, METHODS)
  inm <- METHODS[method]
  if (is.na(method)) 
    stop("invalid distance method")
  if (method == -1) 
    stop("ambiguous distance method")
  if (!method %in% c(1, 2, 6, 16) && any(rowSums(x, na.rm = TRUE) == 
                                         0)) 
    warning("you have empty rows: their dissimilarities may be meaningless in method ", 
            dQuote(inm))
  if (!method %in% c(1, 2, 3, 6, 16) && any(x < 0, na.rm = TRUE)) 
    warning("results may be meaningless because data have negative entries in method ", 
            dQuote(inm))
  if (method == 11 && any(colSums(x) == 0)) 
    warning("data have empty species which influence the results in method ", 
            dQuote(inm))
  if (method == 6) 
    x <- decostand(x, "range", 2, na.rm = TRUE, ...)
  if (method == 16) 
    x <- veganMahatrans(scale(x, scale = FALSE))
  if (binary) 
    x <- decostand(x, "pa")
  N <- nrow(x <- as.matrix(x))
  if (method %in% c(7, 13, 15) && !identical(all.equal(as.integer(x), 
                                                       as.vector(x)), TRUE)) 
    warning("results may be meaningless with non-integer data in method ", 
            dQuote(inm))
  if (method != 17 & method != 18) {
    d <- .C("veg_distance", x = as.double(x), nr = N, nc = ncol(x), 
            d = double(N * (N - 1)/2), diag = as.integer(FALSE), 
            method = as.integer(method), NAOK = na.rm, PACKAGE = "vegan")$d
    if (method == 10) 
      d <- 2 * d/(1 + d)
    d[d < ZAP] <- 0
    if (any(is.na(d))) 
      warning("missing values in results")
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- paste(if (binary) 
      "binary ", METHODS[method], sep = "")
    attr(d, "call") <- match.call()
    class(d) <- "dist"
  }
  if (method == 17) {
    unifrac <- GUniFrac(x1, tree)
    unifrac <- unifrac$unifracs
    d <- as.dist(unifrac[, , 'd_1'])		# Weighted UniFrac
  }
  if (method == 18) {
    unifrac <- GUniFrac(x1, tree)
    unifrac <- unifrac$unifracs
    d <- as.dist(unifrac[, , 'd_UW'])		# Unweighted UniFrac
  }
  d
}
