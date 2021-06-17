# From https://github.com/richelbilderbeek/pirouette_article/issues/59 :
#
# Add worked examples with mutation rates from
# 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8
library(pirouette)
library(beautier)
library(beastier)
library(testthat)
library(ggplot2)

# Constants
example_no <- 24
rng_seed <- 314
crown_age <- 10
mutation_rates <- c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8)
n_phylogenies_per_mutation_rate <- 5
folder_name <- paste0("example_", example_no)
is_testing <- is_on_ci()
if (is_testing) {
  mutation_rates <- c(100, 248)
  n_phylogenies_per_mutation_rate <- 2
}
n_mutation_rates <- length(mutation_rates)
n_pir_params <- n_mutation_rates * n_phylogenies_per_mutation_rate

# Create phylogenies
phylogenies <- list()
for (i in seq_len(n_phylogenies_per_mutation_rate)) {
  set.seed(rng_seed)
  phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)
  phylogenies[[i]] <- phylogeny
}
# 1 2 3 1 2 3
phylogenies <- rep(phylogenies, n_mutation_rates)
expect_equal(n_pir_params, length(phylogenies))

# Create pirouette parameter sets
pir_paramses <- create_std_pir_paramses(
  n = n_pir_params,
  folder_name = folder_name
)
expect_equal(length(pir_paramses), length(phylogenies))
if (is_testing) {
  pir_paramses <- shorten_pir_paramses(pir_paramses)
}

# Set the alignment lengths
# 1 1 1 2 2 2
mutation_rateses <- rep(
  mutation_rates, each = n_phylogenies_per_mutation_rate
)
expect_equal(length(mutation_rateses), length(pir_paramses))
for (i in seq_along(mutation_rateses)) {

  pir_paramses[[i]]$alignment_params$sim_tral_fun <- get_sim_tral_with_std_nsm_fun(
    mutation_rate = mutation_rateses[[i]],
    site_model = beautier::create_jc69_site_model()
  )
  pir_paramses[[i]]$twinning_params$sim_twal_fun <- get_sim_twal_same_n_muts_fun(
    mutation_rate = mutation_rateses[[i]],
    max_n_tries = 1000
  )
}

# Do the runs
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

# Save summary
pir_plots(pir_outs) +
  ggtitle(paste("Number of pir_params: ", n_pir_params)); ggsave("errors.png", width = 7, height = 7)

# Save useful
for (i in seq_along(mutation_rates)) {
  n <- mutation_rates[i]
  from_index <- ((i - 1) * n_phylogenies_per_mutation_rate) + 1
  to_index <- ((i - 1) * n_phylogenies_per_mutation_rate) + n_phylogenies_per_mutation_rate
  pir_plots(
    pir_outs = pir_outs[from_index:to_index]
  ) + ggtitle(
      paste("Mutation rate: ", n, ", number of replicates:", n_phylogenies_per_mutation_rate)
  ); ggsave(filename = paste0("errors_", i, ".png"), width = 7, height = 7)
}

# Save
expect_equal(length(pir_paramses), length(pir_outs))
expect_equal(length(pir_paramses), length(phylogenies))
for (i in seq_along(pir_outs)) {
  pir_save(
    phylogeny = phylogenies[[i]],
    pir_params = pir_paramses[[i]],
    pir_out = pir_outs[[i]],
    folder_name = dirname(pir_paramses[[i]]$alignment_params$fasta_filename)
  )
}
