testthat::test_that(
  desc = "DetermineIsomainPeaks works", 
  code = {
    apci_spectrum <- InterpretMSSpectrum::apci_spectrum
    out1 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=apci_spectrum, ionization="APCI")
    testthat::expect_length(out1, 5)
    testthat::expect_equal(out1, c(202.0893, 246.1339, 274.1288, 348.1477, 364.1789))
    
    esi_spectrum <- InterpretMSSpectrum::esi_spectrum
    out2 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=esi_spectrum, ionization="ESI")
    testthat::expect_length(out2, 5)
    testthat::expect_equal(attr(out2, "pot_mh"), 263.05340)
    
    # exclude double charged peaks in ESI because IMSS can't handle these while findMAIN can
    inp <- data.frame("mz"=c(100,200,200.5,201), "int"=c(100,1000,100,10))
    out3 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="ESI")
    testthat::expect_equal(out3[1], 100)
    out4 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="APCI")
    testthat::expect_equal(out4, c(100,200))
    
    # in ESI peaks below precursor are excluded
    inp <- data.frame("mz"=c(100,200,400), "int"=c(100,100,1000))
    out5 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="ESI")
    testthat::expect_length(out5, 3)
    out6 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="ESI", precursor=200)
    testthat::expect_length(out6, 2)
    
    # limit isomain peak number
    out7 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="ESI", limit=2)
    testthat::expect_length(out7, 2)
    
    # check that at least one peak is returned
    out8 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="ESI", limit=0)
    testthat::expect_length(out8, 1)
    
    # check that at least one peak is returned
    out9 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=apci_spectrum, ionization="APCI", limit=0)
    testthat::expect_length(out9, 1)
    
    out10 <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=apci_spectrum, ionization="APCI", int_cutoff = 0, limit=5)
    testthat::expect_equal(out10, c(202.0893, 246.1339, 274.1288, 348.1477, 364.1789))
    
    # removal of TMS adducts when smaller in int than precursor
    inp <- rbind(apci_spectrum, c(364.1789 + 72.0395, 10^6))
    out <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="APCI", limit=1)
    testthat::expect_equal(out, 364.1789 + 72.0395)
    inp <- rbind(apci_spectrum, c(364.1789 + 72.0395, 10^5))
    out <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="APCI", limit=1)
    testthat::expect_equal(out, 364.1789)
    out <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="APCI", limit=1, precursor = 274)
    testthat::expect_equal(out, 274.1288)
    out <- InterpretMSSpectrum:::DetermineIsomainPeaks(spec=inp, ionization="APCI", limit=1, precursor = 246)
    testthat::expect_equal(out, 246.1339)
    
  }
)
