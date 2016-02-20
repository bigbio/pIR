library(pIR)

context("Isoelectric point estimation: Iterative Method")

test_that("Test of the Iterative Methods.. ", {
    # Test for the salomon pkSet
    seq <- "GLPRKILCAIAKKKGKCKGPLKLVCKC"
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "solomon"), 10.5259)

    # Test for rodwell pKSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "rodwell"), 11.4038)

    #Test for emboss pKSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod =  "emboss"), 10.8014)

    #Test for lehninger pkSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "lehninger"), 10.53)

    #Test for grimsley pkSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "grimsley"), 10.4945)

    #Test for patrickios pKSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "patrickios"), 11.2)

    # Test for DtaSelect pkSet
    expect_equal(pIIterative(sequence = seq, pkSetMethod = "DtaSelect"), 10.0249)
}
)
