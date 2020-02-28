# From https://github.com/richelbilderbeek/pirouette_article/issues/59 :
#
# Add worked examples with mutation rates from
# 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8
suppressMessages(library(pirouette))
suppressMessages(library(ggplot2))

################################################################################
# Constants
################################################################################
is_testing <- is_on_travis()

example_no <- 24

seed_to_mutation_rate <- function(rng_seed) {
  mutation_rates <- c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8)
  mutation_rate <- mutation_rates[rng_seed + 1 - 314]
  if (is.na(mutation_rate)) stop("Invalid seed")
  mutation_rate
}
testthat::expect_equal(seed_to_mutation_rate(314), 0.0125)
testthat::expect_equal(seed_to_mutation_rate(315), 0.025)
testthat::expect_equal(seed_to_mutation_rate(316), 0.05)
testthat::expect_equal(seed_to_mutation_rate(317), 0.1)
testthat::expect_equal(seed_to_mutation_rate(318), 0.2)
testthat::expect_equal(seed_to_mutation_rate(319), 0.4)
testthat::expect_equal(seed_to_mutation_rate(320), 0.8)
testthat::expect_error(seed_to_mutation_rate(321))

for (rng_seed in seq(314, 320)) {

  print(rng_seed)

  folder_name <- paste0("example_", example_no, "_", rng_seed)

  set.seed(rng_seed)
  phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)

  pir_params <- create_std_pir_params(folder_name = folder_name)
  pir_params$alignment_params$sim_tral_fun <- get_sim_tral_with_std_nsm_fun(
    mutation_rate = seed_to_mutation_rate(rng_seed),
    site_model = beautier::create_jc69_site_model()
  )
  pir_params$twinning_params$sim_twal_fun <- get_sim_twal_same_n_muts_fun(
    mutation_rate = seed_to_mutation_rate(rng_seed),
    max_n_tries = 1000
  )

  if (is_testing) {
    pir_params <- shorten_pir_params(pir_params)
  }

  errors <- pir_run(
    phylogeny,
    pir_params = pir_params
  )

  utils::write.csv(
    x = errors,
    file = file.path(folder_name, "errors.csv"),
    row.names = FALSE
  )

  pir_plot(errors) +
    ggsave(file.path(folder_name, "errors.png"), width = 7, height = 7)

  pir_to_pics(
    phylogeny = phylogeny,
    pir_params = pir_params,
    folder = folder_name
  )

  pir_to_tables(
    pir_params = pir_params,
    folder = folder_name
  )
}
