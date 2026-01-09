# Load default settings for CHNOSZ
reset()

info <- "[P|T|E].units() do not accept invalid units"
expect_error(P.units("X"), "units of pressure must be either bar or MPa", info = info)
expect_error(T.units("X"), "units of temperature must be either C or K", info = info)
expect_error(E.units("X"), "units of energy must be either cal or J", info = info)
