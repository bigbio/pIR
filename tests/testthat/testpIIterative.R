library(pIR)

context("Isoelectric point estimation: Iterative Method")

test_that("Test of the Iterative Methods.. ", {
    # Test for the solomon pkSet
    seq <- "GLPRKILCAIAKKKGKCKGPLKLVCKC"
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "solomon"), 10.526)

    # Test for rodwell pKSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "rodwell"), 11.404)

    #Test for emboss pKSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod =  "emboss"), 10.801)

    #Test for lehninger pkSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "lehninger"), 10.530)

    #Test for grimsley pkSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "grimsley"), 10.495)

    #Test for patrickios pKSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "patrickios"), 11.200)

    # Test for DtaSelect pkSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "DtaSelect"), 10.025)
}
)
