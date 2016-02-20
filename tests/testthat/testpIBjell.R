library(pIR)

context("Isoelectric point estimation: Bjellvist Method")

test_that("Test of the Bjellvist Methods.. ", {
    # Test for the salomon pkSet
    seq <- "AGAAPYVQAFDSLLAGPVAE"

    expect_equal(pIBjell(sequence = seq, pkSetMethod = "expasy"), 4.0852)

    # Test for skoog pKSet
    expect_equal(pIBjell(sequence = seq, pkSetMethod = "skoog"), 2.9351)

    expect_equal(pIBjell(sequence = seq, pkSetMethod = "calibrated"), 4.0852)

    # Test for skoog pKSet
    expect_equal(pIBjell(sequence = seq, pkSetMethod = "bjell"), 3.0418)
}
)
