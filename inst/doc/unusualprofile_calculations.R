## ----setup, include = FALSE---------------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "ragg_png",
  out.height = "100%",
  out.width = "100%",
  fig.width = 7.2916667,
  fig.height = 7.2916667,
  message = FALSE,
  warning = FALSE,
  cache = F
)

library(dplyr)
library(stringr)
library(ggnormalviolin)
library(ggplot2)
library(purrr)
library(simstandard)
library(unusualprofile)

## ----more-setup, echo=FALSE---------------------------------------------------
myfont <- "sans"

# font_add_google("Titillium Web", "Titillium Web")
theme_set(theme_minimal(16, base_family = myfont))
update_geom_defaults(geom = "text",
                     new = list(family = myfont))
update_geom_defaults(geom = "label",
                     new = list(family = myfont))
update_geom_defaults(geom = "point",
                     new = list(pch = 16))
#Rounds proportions to significant digits both near 0 and 1
PropRound <- function(p, maxDigits = 10) {
  d <- rep(2, length(p))
  pp <- rep(0, length(p))
  for (i in seq(1, length(p))) {
    if (p[i] > 0.99 || p[i] < 0.01) {
      d[i] <-
        min(maxDigits, 1 - 1 - floor(log10(abs(
          ifelse(p[i] < 0.5, p[i], 1 - p[i])
        ))))
    }
    pp[i] <-
      formatC(round(p[i], digits = d[i]), digits = d[i], flag = "")
    if (round(p[i], digits = maxDigits) == 0) {
      pp[i] <- 0
    }
    if (round(p[i], digits = maxDigits) == 1) {
      pp[i] <- 1
    }
    gsub(" ", "", pp[i])
  }
  
  return(gsub(" ", "", pp))
}

numText <- function(x, digits = 2) {
  str_replace(formatC(x, digits = digits, format = "f"), pattern = "\\-", "−")
}

bmatrix <- function(x, digits = 2) {
  x_formatted <- apply(x, 1, formatC, digits = digits, format = "f")
  x_formatted <- apply(x_formatted, 1, str_remove, pattern = "^0")
  x_formatted[x == 0] <- "0"
  x_formatted[x == 1] <- "1"
  paste0("\\begin{bmatrix}\n",
         paste0(
           apply(
             x_formatted,
             MARGIN = 1,
             FUN = paste0,
             collapse = " & "
           ),
           collapse = "\\\\\n"
         ),
         "\n\\end{bmatrix}")
}



## ----worked-model, echo = FALSE-----------------------------------------------
model <- "X =~ 0.95 * X_1 + 
               0.90 * X_2 + 
               0.85 * X_3 + 
               0.60 * X_4"

## ----worked-data, echo = FALSE------------------------------------------------
# Create data.frame, and add the composite score
d <- data.frame(
  X_1 = 2,
  X_2 = 3,
  X_3 = 1,
  X_4 = 2
) %>%
  simstandard::add_composite_scores(m = model)


## ----standardized-scores, echo=FALSE------------------------------------------
# Standardized observed scores
X <- d[1, -5] %>% t %>% as.vector %>% `names<-`(colnames(d)[-5])


## ----one-dimensional, echo = FALSE, fig.align='center', fig.cap="A simple model with standardized loadings", out.width=300----
knitr::include_graphics("One_dimensional.svg")

## ----example-profile, echo = F, fig.cap="Figure 1. Example profile in a standard multivariate normal distribution."----
tibble(
  Variable = paste0("italic(X)[", 1:4, "]"),
  Score = X,
  vjust = c(1.5,-0.5, 1.5, 1.5)
) %>%
  ggplot(aes(Variable, Score)) +
  geom_normalviolin(aes(mu = 0, sigma = 1), fill = "gray", alpha = .4) +
  geom_line(aes(group = 1), color = "firebrick") +
  geom_point(pch = 16, color = "firebrick", size = 2) +
  geom_text(aes(label = Score, vjust = vjust), color = "firebrick") +
  geom_hline(yintercept = d$X_Composite) +
  scale_x_discrete(
    NULL,
    labels = function(l)
      parse(text = l)
  ) +
  scale_y_continuous() +
  annotate(
    geom = "text",
    x = 3.5,
    y = d$X_Composite,
    label = paste0("Composite = ", formatC(d$X_Composite, 2, format = "f")),
    vjust = -.6,
    size = 4.5
  )

## ----X-vector-----------------------------------------------------------------
X <- c(
  X_1 = 2,
  X_2 = 3,
  X_3 = 1,
  X_4 = 2
)

X

## ----variable-names-----------------------------------------------------------
v_observed <- names(X)
v_latent <- "X"
v_composite <- "X_Composite"
v_names <- c(v_observed, v_composite)

## ----matrix-loadings----------------------------------------------------------
lambda <- c(0.95, 0.90, 0.85, 0.60) %>%
  matrix() %>%
  `rownames<-`(v_observed) %>%
  `colnames<-`(v_latent)

lambda

## ----R-X----------------------------------------------------------------------
# Observed Correlations
R_X <- lambda %*% t(lambda) %>%
  `diag<-`(1)

R_X

## ----implied-corelations, results='asis', echo = FALSE------------------------
cat(paste0("$$R_{X} \\approx ", bmatrix(R_X), "$$"))

## ----weight-matrix------------------------------------------------------------
w <-  cbind(diag(4),
            rep(1, 4)) %>%
  `rownames<-`(v_observed) %>%
  `colnames<-`(v_names)
w

## ----sigma-calculations, include= FALSE---------------------------------------
Sigma <- (t(w) %*% R_X %*% w)
Sigma

## ----Sigma-matrix, results='asis', echo = FALSE-------------------------------
cat(paste0("$$\\Sigma \\approx ", bmatrix(Sigma), "$$"))

## ---- ref.label="sigma-calculations", echo = TRUE-----------------------------
Sigma <- (t(w) %*% R_X %*% w)
Sigma

## ----R-matrix, results='asis', echo = FALSE-----------------------------------
cat(paste0("$$R_{All} \\approx ",
           bmatrix(round(cov2cor(Sigma), 2)),
           "$$"))

## ----matrix-purist, echo = TRUE-----------------------------------------------
D_root_inverted <- Sigma %>%
  diag() %>%
  sqrt() %>%
  diag() %>%
  solve() %>%
  `rownames<-`(v_names) %>%
  `colnames<-`(v_names)

R_all <- D_root_inverted %*% Sigma %*% D_root_inverted

R_all

## ----cor2cov------------------------------------------------------------------
# Convert covariance matrix to correlations 
R_all <- cov2cor(Sigma)

R_all

## ----composite-calculations---------------------------------------------------
# Population means of observed variables
mu_X <- rep_along(X, 0)

# Population standard deviations of observed variables
sd_X <- rep_along(X, 1)

# Covariance Matrix
Sigma_X <- diag(sd_X) %*% R_X %*% diag(sd_X)

# Vector of ones
ones <- rep_along(X, 1)

# Standardized composite score
z_C <- c(ones %*% (X - mu_X) / sqrt(ones %*% (Sigma_X) %*% ones))



## ----expected-X---------------------------------------------------------------
# Predicted value of X, given composite score
X_hat <- sd_X * z_C * R_all[v_observed, v_composite] + mu_X

## -----------------------------------------------------------------------------
d_mc <- c(sqrt(t(X - X_hat) %*% solve(Sigma_X) %*% (X - X_hat)))
d_mc

## ----cm-p---------------------------------------------------------------------
# Number of observed variables
k <- length(v_observed)

# Number of composite variables
j <- length(v_composite)

# Cumulative distribution function
p <- pchisq(d_mc ^ 2, df = k - j)
p

