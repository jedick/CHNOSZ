# Load default settings for CHNOSZ
reset()

# Tests added on 20260110

info <- "Expected results for WP02.auxiliary()"
# Reference values for liquid and vapor density from WP02 Table 13.1
expect_equal(round(WP02.auxiliary("rho.liquid", T = 273.16), 3), 999.793, tolerance = 0.005, scale = 1, info = info)
expect_equal(round(WP02.auxiliary("rho.vapor", T = 273.16), 5), 0.00485, info = info)

# TODO: This is a circular comparison - should find independent reference value somewhere
expect_equal(round(WP02.auxiliary("dP.sigma.dT") * 1e8), 18898, info = info)

info <- "Get expected density for water and steam"
msg_info <- "Get expected message"

# Some distance away from the saturation curve
expect_message(rho_steam <- rho.IAPWS95(T = 373.15, P = 1.0, trace = 1), "steam", info = msg_info)
expect_equal(round(rho_steam), 1, info = info)
expect_message(rho_water <- rho.IAPWS95(T = 373.15, P = 1.1, trace = 1), "water", info = msg_info)
expect_equal(round(rho_water), 958, info = info)

# Closer to the saturation curve
expect_message(rho_vapor <- rho.IAPWS95(T = 373.15, P = 1.0141, trace = 1), "close to saturation", info = msg_info)
expect_equal(round(rho_vapor), 1, info = info)
expect_message(rho_liquid <- rho.IAPWS95(T = 373.15, P = 1.0142, trace = 1), "close to saturation", info = msg_info)
expect_equal(round(rho_liquid), 958, info = info)

info <- "Warning for NA density"
expect_warning(rho.IAPWS95(T = 1e15), "problems finding density", info = info)
