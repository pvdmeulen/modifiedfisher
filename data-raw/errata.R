# ---------------------------------------------------------------------------
# Verify the two typos found in Table 2 of van der Meulen, Raymond &
# van der Meulen (2021) Table 2.
#
# Both checks use only base R: Woolf's closed-form asymptotic CI and
# fisher.test(). Neither relies on the modifiedfisher package, so they are an
# independent confirmation that the printed table is wrong, not the package.

source("data-raw/check_table2_typos.R")
