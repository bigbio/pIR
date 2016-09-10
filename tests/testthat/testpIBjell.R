library(pIR)

context("Isoelectric point estimation: Bjellvist Method")

test_that("Test of the Bjellvist Methods.. ", {
    # Test for the salomon pkSet
    seq <- "AGAAPYVQAFDSLLAGPVAE"

    expect_equal(pIBjell(sequence = seq, pkSetMethod = "expasy"), 3.666)

    # Test for skoog pKSet
    expect_equal(pIBjell(sequence = seq, pkSetMethod = "skoog"), 2.935)

    expect_equal(pIBjell(sequence = seq, pkSetMethod = "calibrated"), 4.085)
}
)
