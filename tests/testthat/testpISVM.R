library(pIR)

context("Isoelectric point estimation: SVM methods")

test_that("Test of the SVM Methods...", {
  
  seq <- "AGAAPYVQAFDSLLAGPVAE"
  
  # Test for "default" model
  expect_equal(round(pISVMpeptide(sequence = seq, model="default"), 6), 3.354431)
  
  # Test for "heller" model
  expect_equal(round(pISVMpeptide(sequence = seq, model="heller"),  6), 4.178931)
  
  # Test for "branca" model
  expect_equal(round(pISVMpeptide(sequence = seq, model="branca"),  6), 4.000313)
  
}
)