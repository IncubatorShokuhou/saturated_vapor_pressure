# -------------------------------------------------------------------------------------------
#    saturated_vapor_pressure
# -------------------------------------------------------------------------------------------
# :DESCRIPTION:
#    Calculate the saturation vapor pressure.
#    Origin implementation in IDL:  http://cires1.colorado.edu/~voemel/vp.pro
#
#    Example: Calculate the vapor pressure over liquid water using the WMO formula

#    ```
#    import random
#    temperature_c = random.uniform(0, 100) - 100.
#    P = vapor_pressure(temperature_c, phase = "liquid", formula = "WMO")
#    ```

#    For temperatures above 0 deg C the vapor pressure over liquid water
#    is calculated.
#
#    The optional parameter `phase` changes the calculation to vapor pressure
#    over liquid water over the entire temperature range.

#    The `formula `to be used can be selected with the appropriate keyword
#
#    The current default formulas are `Hyland and Wexler` for liquid and `Goff Gratch` for ice. (hv20040521)
#

import math
import numpy as np
from typing import Optional


def vapor_pressure_elementwise(
    temperature_c: float, phase: str = "liquid", formula: Optional[str] = None
) -> float:
    """[calculate the saturation vapor pressure, element-wise.]

    Args:
        temperature_c (float): [current temperature [degree C]]
        phase (str, optional): [phase of water surface]. Defaults to "liquid".
        formula (Optional[str], optional): [formula used to calculate saturation pressure]. Defaults to None.

    Raises:
        ValueError: [unrecognized formula name over liquid water]
        ValueError: [unrecognized formula name over ice]

    Returns:
        float: [value of calculated saturation vapor pressure [hPa]]
    """

    temperature_k: float = temperature_c + 273.15  # Most formulas use T in [K]
    temperature_triple_point: float = 273.16  # triple point in K

    # Formulas using [C] use the variable temperature_c

    # Calculate saturation pressure over liquid water ----------------------------

    if phase == "liquid":
        # Default uses Hyland and Wexler over liquid.
        # While this may not be the best formula, it is consistent with what Vaisala uses in their system

        if formula is None:
            formula = "HylandWexler"

        if formula == "MartiMauersberger":
            # Goff Gratch formulation
            # Source : Smithsonian Meteorological Tables, 5th edition, p. 350, 1984
            # From original source: Goff and Gratch (1946), p. 107.

            return vapor_pressure_elementwise(
                temperature_c=temperature_c, phase="liquid", formula="GoffGratch"
            )

        elif formula == "HylandWexler":
            # Source Hyland, R. W. and A. Wexler, Formulations for the Thermodynamic Properties of the saturated Phases of H2O from 173.15K to 473.15K, ASHRAE Trans, 89(2A), 500-519, 1983.
            return (
                math.exp(
                    -0.58002206e4 / temperature_k
                    + 0.13914993e1
                    - 0.48640239e-1 * temperature_k
                    + 0.41764768e-4 * temperature_k ** 2.0
                    - 0.14452093e-7 * temperature_k ** 3.0
                    + 0.65459673e1 * math.log(temperature_k)
                )
                / 100.0
            )

        elif formula == "Hardy":
            # Source Hardy, B., 1998, ITS-90 Formulations for Vapor Pressure, Frostpoint temperature, Dewpoint temperature, and Enhancement Factors in the Range −100 to 100° C.
            # The Proceedings of the Third International Symposium on Humidity & Moisture, London, England
            return (
                math.exp(
                    -2.8365744e3 / temperature_k ** 2
                    - 6.028076559e3 / temperature_k
                    + 1.954263612e1
                    - 2.737830188e-2 * temperature_k
                    + 1.6261698e-5 * temperature_k ** 2
                    + 7.0229056e-10 * temperature_k ** 3
                    - 1.8680009e-13 * temperature_k ** 4
                    + 2.7150305 * math.log(temperature_k)
                )
                / 100.0
            )

        elif formula == "Preining":
            # Source : Vehkamaeki, H., M. Kulmala, I. Napari, K. E. J. Lehtinen, C.Timmreck, M. Noppel, and A. Laaksonen (2002), J. Geophys. Res., 107, doi:10.1029/2002JD002184.
            return (
                math.exp(
                    -7235.424651 / temperature_k
                    + 77.34491296
                    + 5.7113e-3 * temperature_k
                    - 8.2 * math.log(temperature_k)
                )
                / 100.0
            )

        elif formula == "Wexler":
            # Wexler, A., Vapor Pressure Formulation for Water in Range 0 to 100 C. A Revision, Journal of Research of the National Bureau of Standards - A. Physics and Chemistry, September - December 1976, Vol. 80A, Nos.5 and 6, 775-785
            # The line of `T**4` was corrected from '-' to '+' following the original citation. (HV 20140819). The change makes only negligible difference
            return (
                math.exp(
                    -0.29912729e4 * temperature_k ** (-2.0)
                    - 0.60170128e4 * temperature_k ** (-1.0)
                    + 0.1887643854e2 * temperature_k ** 0.0
                    - 0.28354721e-1 * temperature_k ** 1.0
                    + 0.17838301e-4 * temperature_k ** 2.0
                    - 0.84150417e-9 * temperature_k ** 3.0
                    + 0.44412543e-12 * temperature_k ** 4.0
                    + 2.858487 * math.log(temperature_k)
                )
                / 100.0
            )

        elif formula == "GoffGratch":
            # Goff Gratch formulation
            # Source : Smithsonian Meteorological Tables, 5th edition, p. 350, 1984
            # From original source: Goff and Gratch (1946), p. 107.
            temperature_steam_point = 373.16  # steam point temperature in K
            # saturation pressure at steam point temperature, normal atmosphere
            e_water_saturation = 1013.246

            return 10.0 ** (
                -7.90298 * (temperature_steam_point / temperature_k - 1.0)
                + 5.02808 * math.log10(temperature_steam_point / temperature_k)
                - 1.3816e-7
                * (
                    10.0 ** (11.344 * (1.0 - temperature_k / temperature_steam_point))
                    - 1.0
                )
                + 8.1328e-3
                * (
                    10.0 ** (-3.49149 * (temperature_steam_point / temperature_k - 1))
                    - 1.0
                )
                + math.log10(e_water_saturation)
            )

        elif formula == "CIMO":
            # Source: Annex 4B, Guide to Meteorological Instruments and Methods of Observation, WMO Publication No 8, 7th edition, Geneva, 2008. (CIMO Guide)
            return 6.112 * math.exp(17.62 * temperature_c / (243.12 + temperature_c))

        elif formula == "MagnusTetens":
            # Source: Murray, F. W., On the computation of saturation vapor pressure, J. Appl. Meteorol., 6, 203-204, 1967.
            # p_saturation = 10.**(7.5*(temperature_c)/(temperature_c+237.5) + 0.7858)         ; Murray quotes this as the original formula and
            return 6.1078 * math.exp(
                17.269388
                * (temperature_k - temperature_triple_point)
                / (temperature_k - 35.86)
            )  # this as the mathematical equivalent in the form of base e.

        elif formula == "Buck":
            # Bucks vapor pressure formulation based on Tetens formula
            # Source: Buck, A. L., New equations for computing vapor pressure and enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981.
            return 6.1121 * math.exp(17.502 * temperature_c / (240.97 + temperature_c))

        elif formula == "Buck2":
            # Bucks vapor pressure formulation based on Tetens formula
            # Source: Buck Research, Model CR-1A Hygrometer Operating Manual, Sep 2001
            return 6.1121 * math.exp(
                (18.678 - temperature_c / 234.5)
                * temperature_c
                / (257.14 + temperature_c)
            )

        elif formula == "WMO":
            # Intended WMO formulation, originally published by Goff (1957)
            # incorrectly referenced by WMO technical regulations, WMO-NO 49, Vol I, General Meteorological Standards and Recommended Practices, App. A, Corrigendum Aug 2000.
            # and incorrectly referenced by WMO technical regulations, WMO-NO 49, Vol I, General Meteorological Standards and Recommended Practices, App. A, 1988.

            return 10.0 ** (
                10.79574 * (1.0 - temperature_triple_point / temperature_k)
                - 5.02800 * math.log10(temperature_k / temperature_triple_point)
                + 1.50475e-4
                * (
                    1.0
                    - 10.0
                    ** (-8.2969 * (temperature_k / temperature_triple_point - 1.0))
                )
                + 0.42873e-3
                * (
                    10.0
                    ** (+4.76955 * (1.0 - temperature_triple_point / temperature_k))
                    - 1.0
                )
                + 0.78614
            )

        elif formula == "WMO2000":
            # WMO formulation, which is very similar to Goff Gratch
            # Source : WMO technical regulations, WMO-NO 49, Vol I, General Meteorological Standards and Recommended Practices, App. A, Corrigendum Aug 2000.

            return 10.0 ** (
                10.79574 * (1.0 - temperature_triple_point / temperature_k)
                - 5.02800 * math.log10(temperature_k / temperature_triple_point)
                + 1.50475e-4
                * (
                    1.0
                    - 10.0
                    ** (-8.2969 * (temperature_k / temperature_triple_point - 1.0))
                )
                + 0.42873e-3
                * (
                    10.0
                    ** (-4.76955 * (1.0 - temperature_triple_point / temperature_k))
                    - 1.0
                )
                + 0.78614
            )

        elif formula == "Sonntag":
            # Source: Sonntag, D., Advancements in the field of hygrometry, Meteorol. Z., N. F., 3, 51-66, 1994.
            return math.exp(
                -6096.9385 * temperature_k ** (-1.0)
                + 16.635794
                - 2.711193e-2 * temperature_k ** 1.0
                + 1.673952e-5 * temperature_k ** 2.0
                + 2.433502 * math.log(temperature_k)
            )

        elif formula == "Bolton":
            # Source: Bolton, D., The computation of equivalent potential temperature, Monthly Weather Report, 108, 1046-1053, 1980. equation (10)
            return 6.112 * math.exp(17.67 * temperature_c / (temperature_c + 243.5))

        elif formula == "Fukuta":
            # Source: Fukuta, N. and C. M. Gramada, Vapor pressure measurement of supercooled water, J. Atmos. Sci., 60, 1871-1875, 2003.
            # This paper does not give a vapor pressure formulation, but rather a correction over the Smithsonian Tables.
            # Thus calculate the table value first, then use the correction to get to the measured value.
            if temperature_c < -39.0:
                return np.nan

            else:
                pressure_saturation_goff_gratch = vapor_pressure_elementwise(
                    temperature_c=temperature_c, phase="liquid", formula="GoffGratch"
                )

                x = temperature_c + 19
                return pressure_saturation_goff_gratch * (
                    0.9992
                    + 7.113e-4 * x
                    - 1.847e-4 * x ** 2.0
                    + 1.189e-5 * x ** 3.0
                    + 1.130e-7 * x ** 4.0
                    - 1.743e-8 * x ** 5.0
                )

        elif formula == "IAPWS":
            # Source: Wagner W. and A. Pruss (2002), The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use, J. Phys. Chem. Ref. Data, 31(2), 387-535.
            # This is the 'official' formulation from the International Association for the Properties of Water and Steam
            # The valid range of this formulation is 273.16 <= T <= 647.096 K and is based on the ITS90 temperature scale.
            # temperature at the critical point [K]
            temperature_critical_point = 647.096
            pressure_critical_point = (
                22.064 * 1e4
            )  # Vapor pressure at the critical point [hPa]
            nu = 1 - temperature_k / temperature_critical_point
            a1 = -7.85951783
            a2 = 1.84408259
            a3 = -11.7866497
            a4 = 22.6807411
            a5 = -15.9618719
            a6 = 1.80122502
            return pressure_critical_point * math.exp(
                temperature_critical_point
                / temperature_k
                * (
                    a1 * nu
                    + a2 * nu ** 1.5
                    + a3 * nu ** 3.0
                    + a4 * nu ** 3.5
                    + a5 * nu ** 4.0
                    + a6 * nu ** 7.5
                )
            )

        elif formula == "MurphyKoop":
            # Source : Murphy and Koop, Review of the vapour pressure of ice and supercooled water for atmospheric applications, Q. J. R. Meteorol. Soc (2005), 131, pp. 1539-1565.
            return (
                math.exp(
                    54.842763
                    - 6763.22 / temperature_k
                    - 4.210 * math.log(temperature_k)
                    + 0.000367 * temperature_k
                    + math.tanh(0.0415 * (temperature_k - 218.8))
                    * (
                        53.878
                        - 1331.22 / temperature_k
                        - 9.44523 * math.log(temperature_k)
                        + 0.014025 * temperature_k
                    )
                )
                / 100.0
            )

        elif formula == "McIDAS":
            # Source : Unknown, Received from Xin Jin <xjin@ssec.wisc.edu>
            A0 = 0.999996876e0
            A1 = -0.9082695004e-2
            A2 = 0.7873616869e-4
            A3 = -0.6111795727e-6
            A4 = 0.4388418740e-8
            A5 = -0.2988388486e-10
            A6 = 0.2187442495e-12
            A7 = -0.1789232111e-14
            A8 = 0.1111201803e-16
            A9 = -0.3099457145e-19
            B = 0.61078e1
            S = A0 + temperature_c * (
                A1
                + temperature_c
                * (
                    A2
                    + temperature_c
                    * (
                        A3
                        + temperature_c
                        * (
                            A4
                            + temperature_c
                            * (
                                A5
                                + temperature_c
                                * (
                                    A6
                                    + temperature_c
                                    * (A7 + temperature_c * (A8 + temperature_c * A9))
                                )
                            )
                        )
                    )
                )
            )
            return B / S ** 8

        else:
            raise ValueError(
                f"Unknown formula for saturation pressure over liquid water surface: {formula}"
            )

    elif phase == "ice":
        # =============================================================================
        # Calculate saturation pressure over ice -------------------------------------
        if temperature_c <= 0:
            # Independent of the formula used for ice,
            # use Hyland Wexler (water) for temperatures above freezing (see above)
            # Source Hyland, R. W. and A. Wexler, Formulations for the Thermodynamic Properties of the saturated Phases of H2O from 173.15K to 473.15K, ASHRAE Trans, 89(2A), 500-519, 1983.
            return (
                math.exp(
                    -0.58002206e4 / temperature_k
                    + 0.13914993e1
                    - 0.48640239e-1 * temperature_k
                    + 0.41764768e-4 * temperature_k ** 2.0
                    - 0.14452093e-7 * temperature_k ** 3.0
                    + 0.65459673e1 * math.log(temperature_k)
                )
                / 100.0
            )

        else:
            # Default uses Goff Gratch over ice.
            # There is little ambiguity in the ice saturation curve. Goff Gratch is widely used.
            if formula is None:
                formula = "GoffGratch"

            if formula == "WMO2000":
                formula = (
                    "WMO"  # There is no typo issue in the WMO formulations for ice
                )

            if formula == "IAPWS":
                # IAPWS does not provide a vapor pressure formulation over ice; use Goff Gratch instead
                formula = "GoffGratch"

            if formula == "MartiMauersberger":
                # Source : Marti, J. and K Mauersberger, A survey and new measurements of ice vapor pressure at temperatures between 170 and 250 K, GRL 20, 363-366, 1993.
                return 10.0 ** (-2663.5 / temperature_k + 12.537) / 100.0

            elif formula == "HylandWexler":
                # Source Hyland, R. W. and A. Wexler, Formulations for the Thermodynamic Properties of the saturated Phases of H2O from 173.15K to 473.15K, ASHRAE Trans, 89(2A), 500-519, 1983.
                return (
                    math.exp(
                        -0.56745359e4 / temperature_k
                        + 0.63925247e1
                        - 0.96778430e-2 * temperature_k
                        + 0.62215701e-6 * temperature_k ** 2.0
                        + 0.20747825e-8 * temperature_k ** 3.0
                        - 0.94840240e-12 * temperature_k ** 4.0
                        + 0.41635019e1 * math.log(temperature_k)
                    )
                    / 100.0
                )

            elif formula == "Wexler":
                # Wexler, A., Vapor pressure formulation for ice, Journal of Research of the National Bureau of Standards-A. 81A, 5-20, 1977.
                return (
                    math.exp(
                        -0.58653696e4 * temperature_k ** (-1.0)
                        + 0.2224103300e2
                        + 0.13749042e-1 * temperature_k ** 1.0
                        - 0.34031775e-4 * temperature_k ** 2.0
                        + 0.26967687e-7 * temperature_k ** 3.0
                        + 0.6918651 * math.log(temperature_k)
                    )
                    / 100.0
                )

            elif formula == "Hardy":
                # Source Hardy, B., 1998, ITS-90 Formulations for Vapor Pressure, Frostpoint temperature, Dewpoint temperature, and Enhancement Factors in the Range 鈥�100 to +100 掳C, The Proceedings of the Third International Symposium on Humidity & Moisture, London, England
                # These coefficients are updated to ITS90 based on the work by Bob Hardy at Thunder Scientific: http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf
                # The difference to the older ITS68 coefficients used by Wexler is academic.
                return (
                    math.exp(
                        -0.58666426e4 * temperature_k ** (-1.0)
                        + 0.2232870244e2
                        + 0.139387003e-1 * temperature_k ** 1.0
                        - 0.34262402e-4 * temperature_k ** 2.0
                        + 0.27040955e-7 * temperature_k ** 3.0
                        + 0.67063522e-1 * math.log(temperature_k)
                    )
                    / 100.0
                )

            elif formula == "GoffGratch":
                # Source : Smithsonian Meteorological Tables, 5th edition, p. 350, 1984
                e_ice_0 = 6.1071  # hPa

                return 10.0 ** (
                    -9.09718 * (temperature_triple_point / temperature_k - 1.0)
                    - 3.56654 * math.log10(temperature_triple_point / temperature_k)
                    + 0.876793 * (1.0 - temperature_k / temperature_triple_point)
                    + math.log10(e_ice_0)
                )

            elif formula == "MagnusTetens":
                # Source: Murray, F. W., On the computation of saturation vapor pressure, J. Appl. Meteorol., 6, 203-204, 1967.
                # p_saturation = 10.**(9.5 * temperature_c/(265.5+temperature_c) + 0.7858)         ; Murray quotes this as the original formula and
                return 6.1078 * math.exp(
                    21.8745584
                    * (temperature_k - temperature_triple_point)
                    / (temperature_k - 7.66)
                )  # this as the mathematical equivalent in the form of base e.

            elif formula == "Buck":
                # Bucks vapor pressure formulation based on Tetens formula
                # Source: Buck, A. L., New equations for computing vapor pressure and enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981.
                return 6.1115 * math.exp(
                    22.452 * temperature_c / (272.55 + temperature_c)
                )

            elif formula == "Buck2":
                # Bucks vapor pressure formulation based on Tetens formula
                # Source: Buck Research, Model CR-1A Hygrometer Operating Manual, Sep 2001
                return 6.1115 * math.exp(
                    (23.036 - temperature_c / 333.7)
                    * temperature_c
                    / (279.82 + temperature_c)
                )

            elif formula == "CIMO":
                # Source: Annex 4B, Guide to Meteorological Instruments and Methods of Observation, WMO Publication No 8, 7th edition, Geneva, 2008. (CIMO Guide)
                return 6.112 * math.exp(
                    22.46 * temperature_c / (272.62 + temperature_c)
                )

            elif formula == "WMO":
                # WMO formulation, which is very similar to Goff Gratch
                # Source : WMO technical regulations, WMO-NO 49, Vol I, General Meteorological Standards and Recommended Practices, Aug 2000, App. A.

                return 10.0 ** (
                    -9.09685 * (temperature_triple_point / temperature_k - 1.0)
                    - 3.56654 * math.log10(temperature_triple_point / temperature_k)
                    + 0.87682 * (1.0 - temperature_k / temperature_triple_point)
                    + 0.78614
                )

            elif formula == "Sonntag":
                # Source: Sonntag, D., Advancements in the field of hygrometry, Meteorol. Z., N. F., 3, 51-66, 1994.
                return math.exp(
                    -6024.5282 * temperature_k ** (-1.0)
                    + 24.721994
                    + 1.0613868e-2 * temperature_k ** 1.0
                    - 1.3198825e-5 * temperature_k ** 2.0
                    - 0.49382577 * math.log(temperature_k)
                )

            elif formula == "MurphyKoop":
                # Source : Murphy and Koop, Review of the vapour pressure of ice and supercooled water for atmospheric applications, Q. J. R. Meteorol. Soc (2005), 131, pp. 1539-1565.
                return (
                    math.exp(
                        9.550426
                        - 5723.265 / temperature_k
                        + 3.53068 * math.log(temperature_k)
                        - 0.00728332 * temperature_k
                    )
                    / 100.0
                )

            elif formula == "McIDAS":
                # Source : Unknown, Received from Xin Jin <xjin@ssec.wisc.edu>
                A0 = 0.7859063157e0
                A1 = 0.3579242320e-1
                A2 = -0.1292820828e-3
                A3 = 0.5937519208e-6
                A4 = 0.4482949133e-9
                A5 = 0.2176664827e-10
                E = A0 + temperature_c * (
                    A1
                    + temperature_c
                    * (
                        A2
                        + temperature_c
                        * (A3 + temperature_c * (A4 + temperature_c * A5))
                    )
                )
                return 10.0 ** E

            else:
                raise ValueError(
                    f"Unknown formula for saturation pressure over ice surface: {formula}"
                )

    else:
        raise ValueError("Phase of water must be either `liquid` or `ice`.")

def vapor_pressure(
    temperature_c: np.ndarray, phase: str = "liquid", formula: Optional[str] = None
) -> np.ndarray:
    """[calculate the saturation vapor pressure.]

    Args:
        temperature_c (float): [current temperature [degree C]]
        phase (str, optional): [phase of water surface]. Defaults to "liquid".
        formula (Optional[str], optional): [formula used to calculate saturation pressure]. Defaults to None.

    Raises:
        ValueError: [unrecognized formula name over liquid water]
        ValueError: [unrecognized formula name over ice]

    Returns:
        np.ndarray: [value of calculated saturation vapor pressure [hPa]]
    """
    if phase == "liquid":
        if formula == "MartiMauersberger":
            print(
                "Marti and Mauersberger don't have a vapor pressure curve over liquid. Using Goff Gratch instead"
            )
            formula = "GoffGratch"
        return np.vectorize(pyfunc=vapor_pressure_elementwise)(
            temperature_c=temperature_c, phase="liquid", formula=formula
        )

    elif phase == "ice":
        if formula == "IAPWS":
            print(
                "IAPWS does not provide a vapor pressure formulation over ice; use Goff Gratch instead"
            )
            formula = "GoffGratch"

        return np.vectorize(pyfunc=vapor_pressure_elementwise)(
            temperature_c=temperature_c, phase="ice", formula=formula
        )

    else:
        raise ValueError("Phase of water must be either `liquid` or `ice`.")