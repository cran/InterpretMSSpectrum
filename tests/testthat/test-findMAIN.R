testthat::test_that(
  desc = "findMAIN works", 
  code = {
    esi_spectrum <- InterpretMSSpectrum::esi_spectrum
    fmr <- InterpretMSSpectrum::findMAIN(esi_spectrum)
    
    # returns a object of class 'findMAIN'
    testthat::expect_true(inherits(fmr, "findMAIN"))
    
    # result should have 96 hypothesis tested (and be of length 96 therefore)
    testthat::expect_length(fmr, 96)
    
    # compute summary of fmr
    fmr_sum <- summary(fmr)
    
    # summary has 10 columns
    testthat::expect_equal(ncol(fmr_sum), 10)
    
    # summary has results for 3 main adduct hypotheses
    testthat::expect_true(all(fmr_sum[,"adducthyp"] %in% c("[M+H]+","[M+Na]+","[M+K]+")))
    
    # compute an alternative for only 4 peaks out of spec using only one adduct hypothesis
    fmr_sum2 <- summary(InterpretMSSpectrum::findMAIN(esi_spectrum[6:9,], adducthyp = "[M+H]+"))
    
    # summary2 has 2 rows
    testthat::expect_equal(nrow(fmr_sum2), 2)
    
    # summary2 has only M+H adduct hypotheses
    testthat::expect_true(all(fmr_sum2[,"adducthyp"] == c("[M+H]+")))
    
    # set up a spectrum containing a double charged peak
    spec <- data.frame(mz = c(372.1894, 372.6907, 373.1931, 380.1234), int = c(100, 40, 8, 2))
    
    # allow a double charged adduct hypothesis (not standard)
    fmr <- InterpretMSSpectrum::findMAIN(spec, adducthyp = c("[M+H]+", "[M+2H]2+"))
    fmr_sum3 <- summary(fmr)
    
    # check that correct adduct hyp is ranked first
    testthat::expect_true(fmr_sum3[1,2]=="[M+2H]2+")
    
    # add the correct M+H to this spectrum as a minor peak
    spec <- rbind(spec, c(742.3648+1.007, 10))
    fmr_sum4 <- summary(InterpretMSSpectrum::findMAIN(spec, adducthyp = c("[M+H]+", "[M+2H]2+")))
    
    # summary4 has 5 rows because one redundant hypothesis is removed
    testthat::expect_equal(nrow(fmr_sum4), 5)
    
    # compare specific hypotheses manually
    # get correct result and check total_score
    fmr_sum5 <- summary(InterpretMSSpectrum::findMAIN(spec, adductmz = 743.3718, adducthyp = "[M+H]+"))
    testthat::expect_equal(fmr_sum5[,10], 0.95)
    # enforce wrong result and check total_score
    fmr_sum6 <- summary(InterpretMSSpectrum::findMAIN(spec, adductmz = 743.3718, adducthyp = "[M+2H]2+"))
    testthat::expect_equal(fmr_sum6[,10], 0.51)
  }
)

testthat::test_that(
  desc = "findMAIN plot returns expected result",
  code = {
    # avoid creating a Rplots.pdf in testthat folder
    pdf(NULL)
    esi_spectrum <- InterpretMSSpectrum::esi_spectrum
    fmr <- InterpretMSSpectrum::findMAIN(esi_spectrum)
    vdiffr::expect_doppelganger(
      title = "findMAIN_Plot_01",
      fig = function() plot(fmr)
    )
  }
)

testthat::test_that(
  desc = "findMAIN print returns expected result",
  code = {
    esi_spectrum <- InterpretMSSpectrum::esi_spectrum
    fmr <- InterpretMSSpectrum::findMAIN(esi_spectrum)
    suppressMessages({
      testthat::expect_output(print(fmr), "charge")
    })
  }
)

testthat::test_that(
  desc = "findMAIN works for negative charged data", 
  code = {
    
    ads <- c("[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-", 
             "[2M-H]-", "[2M+Na-2H]-", "[2M+K-2H]-", "[2M+CH2O2- H]-",
             "[3M-H]-", "[3M+Na-2H]-", "[3M+K-2H]-", "[3M+CH2O2- H]-")
    
    spec <- structure(list(mz = c(59.0138, 71.0137, 85.0292, 89.0244, 101.0241, 
                                  107.0354, 113.0242, 114.0275, 119.0351, 120.0387, 143.0349, 149.0456, 
                                  161.0453, 162.049, 175.04, 179.0557, 180.0593, 181.0605, 225.0615, 
                                  226.0651, 227.0661, 232.9754, 247.043, 292.9965, 303.9916, 309.0141, 
                                  310.093, 311.0968, 312.0984, 321.0804, 337.0548, 356.0984, 357.1019, 
                                  358.1041, 373.0879, 374.0912, 378.0792, 381.1013, 397.0754, 422.0384, 
                                  423.0268, 424.0337, 424.0793, 425.0368, 425.0797, 440.0517, 484.0545, 
                                  512.1381, 513.1416, 543.0684, 544.0746, 558.1436, 606.5348, 607.0366, 
                                  643.1745, 644.1766, 645.1514, 646.1548, 657.17, 659.1506, 660.1516, 
                                  665.1578, 666.163, 674.1047, 675.1103, 676.1144, 689.174, 691.1566, 
                                  692.16, 719.1236, 720.1125, 721.1175, 722.1206, 723.1237, 733.2031, 
                                  733.7053, 736.1803, 738.2728, 798.7207, 799.2218, 799.7195, 814.6903, 
                                  815.1911, 821.723, 822.2251, 822.7256, 836.7164, 954.2711, 954.7715, 
                                  955.2724, 955.7684, 956.2506, 957.2531, 958.2555, 962.2551, 962.7569, 
                                  970.2402, 970.7412, 974.2605, 976.2541, 977.2636, 977.7725, 985.2043, 
                                  986.2096, 987.2128, 1007.187, 1109.8184, 1110.3202, 1110.8202
    ), int = c(2103.63924946167, 3562.9507269249, 11762.0336577921, 
               27681.0835028321, 8688.67599721221, 2538.51227944622, 38220.8995936023, 
               2061.80244027094, 40795.7856726914, 1653.01052157474, 12317.3976616045, 
               6480.27739293513, 32869.3698615668, 2113.09337687844, 1183.08038920853, 
               989480.409255996, 66252.3677592906, 12622.9110591638, 171269.593577525, 
               12521.4746941321, 2909.9100904207, 142.726446390455, 31304.3949157096, 
               80.4776723752696, 85.0640673818838, 2363.69089814068, 212566.41069358, 
               32620.6006573411, 5271.19705275049, 122.14521094642, 2630.77948106963, 
               285427.686649552, 47096.0815937793, 8164.5868577798, 128951.658610142, 
               19245.3177148242, 18951.7499754403, 489.474571988212, 4524.78620358222, 
               230.264803027306, 216.039594767078, 289.95029998503, 47564.9243008665, 
               1275.2437293667, 11234.5770149179, 23447.1150796509, 456.772806162486, 
               1308.30693959136, 5111.64488634635, 4170.5007412817, 3428.76909968807, 
               1621.97682432071, 3242.31398225046, 1930.28587896828, 3449.65139632296, 
               1299.62373770681, 15734.4584055888, 4929.39290109019, 3345.05755927916, 
               20861.8830244032, 7067.96224581355, 7127.38572799359, 3045.93390602148, 
               8007.76672924208, 13293.4360007774, 4325.09312163993, 3300.23745819813, 
               22820.9749628119, 7676.01459869228, 4841.80502761374, 9597.72173603877, 
               79322.314205121, 25406.0700572985, 6906.98581612917, 476.241912973919, 
               362.524687348738, 1514.55469310517, 2365.33760578423, 2427.62247542958, 
               2260.02106023082, 2037.80313303156, 2640.74126462242, 2179.85491619748, 
               3267.42390941161, 2836.70048525699, 1549.1351299129, 1852.05763080516, 
               7254.39542583341, 4924.4679938275, 4704.68062595821, 2062.33543451824, 
               21950.7409323004, 10432.9446180286, 4042.32978225466, 1710.69669128373, 
               1394.75255358217, 3237.3597910819, 2652.15613816062, 2985.1276013794, 
               5539.98500001492, 5008.33165132181, 2751.00677937685, 15390.0325516664, 
               18784.1184381019, 8095.68377920769, 3252.01168141609, 2318.97801620158, 
               2662.15861453051, 2137.26181104054)), row.names = c(40L, 46L, 
                                                                   3L, 64L, 24L, 73L, 53L, 22L, 32L, 25L, 83L, 56L, 100L, 104L, 
                                                                   76L, 27L, 47L, 63L, 6L, 54L, 77L, 65L, 105L, 51L, 41L, 107L, 
                                                                   5L, 8L, 9L, 34L, 4L, 48L, 74L, 88L, 80L, 10L, 106L, 15L, 1L, 
                                                                   66L, 33L, 14L, 109L, 7L, 108L, 103L, 67L, 60L, 2L, 97L, 92L, 
                                                                   68L, 57L, 94L, 28L, 42L, 58L, 59L, 30L, 93L, 69L, 43L, 21L, 98L, 
                                                                   61L, 62L, 70L, 101L, 96L, 91L, 102L, 84L, 72L, 95L, 12L, 71L, 
                                                                   52L, 79L, 75L, 31L, 44L, 81L, 89L, 99L, 85L, 26L, 86L, 45L, 36L, 
                                                                   19L, 23L, 18L, 13L, 11L, 20L, 17L, 82L, 35L, 78L, 55L, 16L, 90L, 
                                                                   37L, 29L, 38L, 87L, 49L, 50L, 39L), class = "data.frame")
    fmr <- InterpretMSSpectrum::findMAIN(
      spec, 
      rules = ads, 
      adducthyp = ads[grep("[M", ads, fixed = TRUE)], 
      ionmode = 'negative', 
      mzabs = 0.005, 
      ppm = 10,
      mainpkthr = 0.15
    )
    
    # compute summary of fmr
    fmr_sum <- summary(fmr)
    
    # summary has 10 columns
    testthat::expect_equal(ncol(fmr_sum), 10)
    
    # check that correct adduct hyp is ranked first and the total_score is consistent
    testthat::expect_true(fmr_sum[1,"adducthyp"]=="[M-H]-")
    testthat::expect_true(fmr_sum[1,"total_score"]==0.58)
    
    # check that all neutral masses are positive
    testthat::expect_true(all(fmr_sum[,"neutral_mass"]>0))
    
    # check that supp_isos is finite
    testthat::expect_true(all(is.finite(fmr_sum[,"supp_isos"])))
    
    # check that medppm is finite when more than 1 adduct is found
    testthat::expect_true(all(is.finite(fmr_sum[fmr_sum[,"adducts_explained"]>=2,"medppm"])))
    
  }
)