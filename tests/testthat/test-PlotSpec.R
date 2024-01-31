testthat::test_that(
  desc = "PlotSpec plot returns expected result for APCI",
  code = {
    inp <- InterpretMSSpectrum::apci_spectrum
    # avoid creating a Rplots.pdf in testthat folder
    pdf(NULL)
    vdiffr::expect_doppelganger(
      title = "PlotSpec_01",
      fig = function() InterpretMSSpectrum::PlotSpec(x = inp, ionization="APCI")
    )
  }
)

testthat::test_that(
  desc = "PlotSpec plot returns expected result for ESI",
  code = {
    inp <- InterpretMSSpectrum::esi_spectrum
    substitutions <- InterpretMSSpectrum::param.default$substitutions
    neutral_losses <- data.frame("Name"="Test", "Formula"="C_X", "Mass"=291.049-263.0534)
    # avoid creating a Rplots.pdf in testthat folder
    pdf(NULL)
    vdiffr::expect_doppelganger(
      title = "PlotSpec_02",
      fig = function() InterpretMSSpectrum::PlotSpec(x = inp, ionization="ESI", xlim = c(260, 300), substitutions = substitutions, neutral_losses = neutral_losses, precursor = 291)
    )
  }
)
