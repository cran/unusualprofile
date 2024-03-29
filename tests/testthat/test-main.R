# Make data ----

library(testthat)
library(simstandard)
library(dplyr)
library(unusualprofile)
# lavaan model with three indicators of a latent variable
model <- "
X =~ 0.7 * X_1 + 0.5 * X_2 + 0.8 * X_3
Y =~ 0.7 * Y_1 + 0.5 * Y_2 + 0.8 * Y_3
Y ~ 0.6 * X
"

# Randomly generated case
d <- sim_standardized(
  model,
  n = 1,
  observed = TRUE,
  latent = TRUE,
  errors = FALSE,
  composites = TRUE
)

# Model-implied correlation matrix
v_dep <- c("X_1", "X_2", "X_3",
           "Y_1", "Y_2", "Y_3")
v_ind_composites <- c("X_Composite", "Y_Composite")
v_ind <- NULL
v_all <- c(v_ind, v_ind_composites, v_dep)
R <- sim_standardized_matrices(model)$Correlations$R_all[v_all,
                                                         v_all,
                                                         drop = FALSE]


cond_maha(
  data = d,
  R = R,
  v_dep = v_dep,
  v_ind_composites = v_ind_composites
)

# Works  and Input ----
test_that("Example works", {
  expect_silent({
    dcm <- cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )
    dcm
    plot(dcm)
  })

})



test_that("data as matrix okay", {
  expect_silent(cond_maha(
    data = as.matrix(d),
    R = R,
    v_dep = v_dep,
    v_ind_composites = v_ind_composites
  ))
})

test_that("data as vector okay", {
  expect_silent(cond_maha(
    as.matrix(d)[1,],
    R = R,
    v_dep = v_dep,
    v_ind_composites = v_ind_composites
  ))
})

test_that("R is a correlation matrix", {
  R1 <- R
  R1[1, 2] <- 1.1
  R1[2, 1] <- 1.1
  expect_error(
    cond_maha(
      data = d,
      R = R1,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    ),
    "Some values of R are outside the range of 1 and -1."
  )
  R1 <- R
  R1[1, 2] <- -1.1
  R1[2, 1] <- -1.1
  expect_error(
    cond_maha(
      data = d,
      R = R1,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    ),
    "Some values of R are outside the range of 1 and -1."
  )
  R1 <- R
  R1[1, 1] <- 0.5
  expect_error(
    cond_maha(
      data = d,
      R = R1,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    ),
    "R has values on its diagonal that are not ones."
  )
  R1 <- R
  R1[1, 2] <- 0.5
  R1[2, 1] <- 0.6
  expect_error(
    cond_maha(
      data = d,
      R = R1,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    ),
    "R is not symmetric"
  )
})

# Variable names ----

test_that("v_ind, v_ind_composites, and v_dep are in colnames of data", {
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind = "fred",
      v_ind_composites = v_ind_composites
    ),
    "Some variables in v_ind are not in data"
  )
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = c(v_dep, "fred"),
      v_ind_composites = v_ind_composites
    ),
    "Some variables in v_dep are not in data"
  )
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = c(v_ind_composites, "fred")
    ),
    "Some variables in v_ind_composites are not in data"
  )
})


test_that("v_ind, v_ind_composites, and v_dep are in colnames of R", {
  d1 <- d %>% mutate(fred = 1)

  expect_error(
    cond_maha(
      data = d1,
      R = R,
      v_dep = v_dep,
      v_ind = "fred",
      v_ind_composites = v_ind_composites
    ),
    "Some variables in v_ind are not in R: fred"
  )

  expect_error(
    cond_maha(
      data = d1,
      R = R,
      v_dep = c(v_dep, "fred"),
      v_ind = NULL,
      v_ind_composites = v_ind_composites
    ),
    "Some variables in v_dep are not in R: fred"
  )

  # v_ind_composites are null
  expect_no_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = c("Y_1", "Y_2"),
      v_ind = c("X_1", "X_2"),
      v_ind_composites = NULL
    )
  )

  expect_error(
    cond_maha(
      data = d1,
      R = R,
      v_dep = v_dep,
      v_ind_composites = c(v_ind_composites, "fred")
    ),
    "Some variables in v_ind_composites are not in R: fred"
  )
})

test_that("v_ind and v_dep have no overlapping names", {
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind = "X_1",
      v_ind_composites = v_ind_composites
    )
  )
  expect_error(cond_maha(
    data = d,
    R = R,
    v_dep = c(v_dep, "X"),
    v_ind_composites = v_ind_composites
  ))
  expect_error(cond_maha(
    data = d,
    R = R,
    v_dep = v_dep,
    v_ind_composites = c(v_ind_composites, "X_1")
  ))
})

test_that("Sample stats need at least 3 rows of data.", {
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      use_sample_stats = TRUE
    )
  )

  expect_no_error(
    cond_maha(
      data = sim_standardized(
        model,
        n = 100,
        observed = TRUE,
        latent = TRUE,
        errors = FALSE,
        composites = TRUE
      ),
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      use_sample_stats = TRUE
    )
  )


  expect_error(plot(cond_maha(
    data = sim_standardized(
      model,
      n = 2,
      observed = TRUE,
      latent = TRUE,
      errors = FALSE,
      composites = TRUE
    ),
    R = R,
    v_dep = v_dep,
    v_ind_composites = v_ind_composites
  )), "Can only plot one case at a time")

  expect_error(plot(cond_maha(
    data = sim_standardized(
      model,
      n = 2,
      observed = TRUE,
      latent = TRUE,
      errors = FALSE,
      composites = TRUE
    ),
    R = R,
    v_dep = v_dep
  )), "Can only plot one case at a time")

})

# Singularity ----

test_that("Measures not singular", {
  R_singular_dep <- R
  R_singular_dep["Y_1", "Y_2"] <- R_singular_dep["Y_2", "Y_1"] <- 1
  R_singular_ind <- R
  R_singular_ind["X_Composite", "Y_Composite"] <-
    R_singular_ind["Y_Composite", "X_Composite"] <- 1
  expect_error(
    cond_maha(
      data = d,
      R = R_singular_dep,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )
  )
  expect_error(
    cond_maha(
      data = d,
      R = R_singular_ind,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )
  )

  expect_no_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = c("Y_1", "Y_2", "Y_3"),
      v_ind = "X_Composite",
      v_ind_composites = "Y_Composite"
    )
  )

  expect_true(is_singular(matrix(1, 2, 2)))
})

# Variable Metrics ----

test_that("Metric does not matter", {
  # All have different metric
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )$dCM,
    cond_maha(
      data = d %>% mutate_all(function(x) {
        x * 15 + 100
        }
        ),
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = 100,
      sigma = 15
    )$dCM
  )
  # 1 dep has a different metric
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )$dCM,
    cond_maha(
      data = d %>% mutate(X_1 = X_1 * 15 + 100),
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = c(X_1 = 100, X_2 = 0, X_3 = 0, Y_1 = 0, Y_2 = 0, Y_3 = 0,
             X_Composite = 0, Y_Composite = 0),
      sigma = c(X_1 = 15, X_2 = 1, X_3 = 1, Y_1 = 1, Y_2 = 1, Y_3 = 1,
                X_Composite = 1, Y_Composite = 1)
    )$dCM
  )

  # 1 composite has a different metric
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )$dCM,
    cond_maha(
      data = d %>% mutate(X_Composite = X_Composite * 15 + 100),
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = c(X_1 = 0, X_2 = 0, X_3 = 0, Y_1 = 0, Y_2 = 0, Y_3 = 0,
             X_Composite = 100, Y_Composite = 0),
      sigma = c(X_1 = 1, X_2 = 1, X_3 = 1, Y_1 = 1, Y_2 = 1, Y_3 = 1,
                X_Composite = 15, Y_Composite = 1)
    )$dCM
  )

  # 1 ind has a different metric
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = c("Y_1", "Y_2", "Y_3"),
      v_ind = c("X_1", "X_2", "X_3"),
    )$dCM,
    cond_maha(
      data = d %>% mutate(X_3 = X_3 * 15 + 100),
      R = R,
      v_dep = c("Y_1", "Y_2", "Y_3"),
      v_ind = c("X_1", "X_2", "X_3"),
      mu = c(X_1 = 0, X_2 = 0, X_3 = 100, Y_1 = 0, Y_2 = 0, Y_3 = 0),
      sigma = c(X_1 = 1, X_2 = 1, X_3 = 15, Y_1 = 1, Y_2 = 1, Y_3 = 1)
    )$dCM
  )

})

# Rounding ----
test_that("Propround", {
  expect_equal(
    proportion_round(
      c(
        -1,
        -.0001,
        0,
        .000012,
        0.012,
        0.12,
        0.98,
        0.99,
        0.991,
        0.999,
        0.9991,
        1,
        1.001,
        2,
        NA
      )
    ),
    c(
      "-1.00",
      "-0.00",
      ".00",
      ".000012",
      ".01",
      ".12",
      ".98",
      ".99",
      ".991",
      ".999",
      ".9991",
      "1.00",
      "1.00",
      "2.00",
      NA_character_
    )
  )

  expect_equal(
    proportion_round(
      c(
        -1,-.00012,
        0,
        .00001234,
        0.01234,
        0.1234,
        0.9111,
        0.99,
        0.991111,
        0.999,
        0.9991111,
        1,
        1.001234,
        2,
        NA
      ),
      digits = 3
    ),
    c(
      "-1.000",
      "-0.000",
      ".000",
      ".000012",
      ".012",
      ".123",
      ".911",
      ".990",
      ".991",
      ".999",
      ".9991",
      "1.000",
      "1.001",
      "2.000",
      NA_character_
    )
  )

  expect_equal(proportion2percentile(c(0.001, 0.01, .5, .99, .992)),
               c(".1", "1", "50", "99", "99.2"))

  expect_equal(
    proportion2percentile(c(0.001, 0.01, .5, .99, .992),
                          add_percent_character = TRUE),
    c(".1%", "1%", "50%", "99%", "99.2%")
  )
})

test_that("Labeling", {
  expect_equal(
    p2label(pnorm(
      seq(60, 140, 10),
      mean = 100,
      sd = 15
    )),
    c(
      "Extremely Low",
      "Very Low",
      "Low",
      "Low Average",
      "Average",
      "High Average",
      "High",
      "Very High",
      "Extremely High"
    )
  )
})

# Independent and Independent Composite Together ----

test_that("ind and ind_composite together", {
  expect_silent(cond_maha(
    data = d,
    R = R,
    v_dep = c("Y_1", "Y_2", "Y_3"),
    v_ind = c("X_1", "X_2", "X_3"),
    v_ind_composites = "Y_Composite"
  ))
})

# Mu and sigma names not assigned

test_that("mu and sigma names not assigned", {
  # mu and sigma are null
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )$dCM,
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = NULL,
      sigma = NULL
    )$dCM
  )



  # mu and sigma are scalars
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )$dCM,
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = 0,
      sigma = 1
    )$dCM
  )

  # mu and sigma are unnamed vectors with equal values
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )$dCM,
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = rep(0, 8),
      sigma = rep(1, 8)
    )$dCM
  )

  # mu and sigma are unnamed vectors and not the same
  expect_equal(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites
    )$dCM,
    cond_maha(
      data = d %>% mutate(X_1 = X_1 * 15 + 100),
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = c(100, rep(0, 7)),
      sigma = c(15, rep(1, 7))
    )$dCM
  )

  # Too many mus
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      mu = rep(0, 10)
    )$dCM,
    regexp = paste0(
      "There are 10 means in mu, ",
      "but 8 columns in the data. ",
      "Supply 1 mean for each variable."
    )
  )

  # Too many sigmas
  expect_error(
    cond_maha(
      data = d,
      R = R,
      v_dep = v_dep,
      v_ind_composites = v_ind_composites,
      sigma = rep(0,10)
    )$dCM,regexp = paste0("There are 10 means in sigma, ",
                          "but 8 columns in the data. ",
                          "Supply 1 mean for each variable.")
  )



})

# Mahalanobis
test_that("No predictors", {
  expect_no_error({
    dm <- cond_maha(data = d,
                    R = R,
                    v_dep = v_dep)
    dm
    plot(dm)
  })
})




