import math


def saturated_vapor_pressure(T: float, surface="water", method="Goff-Grattch") -> float:
    """[calculates the saturated vapor pressure]
    See:
        http://cires1.colorado.edu/~voemel/vp.html

    Args:
        T (float): [temperature [K]]
        surface (str, optional): [surface type (water or ice)]. Defaults to "water".
        method (str, optional): [methods of calculation]. Defaults to "Goff-Grattch".

    Returns:
        float: [saturated vapor pressure [hPa]]
    """
    if surface == "water":
        if (
            method == "Goff-Gratch"
        ):  # (Smithsonian Tables, 1984, after Goff and Gratch, 1946)
            Log10ew = (
                -7.90298 * (373.16 / T - 1)
                + 5.02808 * math.log10(373.16 / T)
                - 1.3816 * 10 ** -7 * (10 ** (11.344 * (1 - T / 373.16)) - 1)
                + 8.1328 * 10 ** -3 * (10 ** (-3.49149 * (373.16 / T - 1)) - 1)
                + math.log10(1013.246)
            )
            return 10 ** Log10ew

        elif (
            method == "CIMO"
        ):  # Guide to Meteorological Instruments and Methods of Observation (CIMO Guide)
            tc = T - 273.15
            return 6.112 * math.exp(17.67 * tc / (243.12 + tc))

        elif method == "WMO2012":  # Goff, J. A. Saturation pressure of water on the new Kelvin temperature scale, Transactions of the American society of heating and ventilating engineers, pp 347-354, presented at the semi-annual meeting of the American society of heating and ventilating engineers, Murray Bay, Que. Canada, 1957.; World Meteorological Organization, Technical Regulations, Basic Documents No. 2, Volume I - General meteorological standards and recommended practices, Appendix A, WMO-No. 49, Geneva 2011, updated 2012.
            Log10ew = (
                10.79574 * (1 - 273.16 / T)
                - 5.02800 * math.log10(T / 273.16)
                + 1.50475 * 10 ** -4 * (1 - 10 ** (-8.2969 * (T / 273.16 - 1)))
                + 0.42873 * 10 ** -3 * (10 ** (4.76955 * (1 - 273.16 / T)) - 1)
                + 0.78614
            )
            return 10 ** Log10ew

        # Hyland, R. W. and A. Wexler, Formulations for the Thermodynamic Properties of the saturated Phases of H2O from 173.15K to 473.15K, ASHRAE Trans, 89(2A), 500-519, 1983.
        elif method == "Hyland-Wexler":
            Logew = (
                -0.58002206 * 10 ** 4 / T
                + 0.13914993 * 10 ** 1
                - 0.48640239 * 10 ** -1 * T
                + 0.41764768 * 10 ** -4 * T ** 2
                - 0.14452093 * 10 ** -7 * T ** 3
                + 0.65459673 * 10 ** 1 * math.log(T)
            )
            return math.exp(Logew) / 100

        elif method == "Hardy":  # Hardy, B., 1998, ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature, and Enhancement Factors in the Range –100 to +100 °C, The Proceedings of the Third International Symposium on Humidity & Moisture, London, England
            Logew = (
                -2.8365744 * 10 ** 3 / T ** 2
                - 6.028076559 * 10 ** 3 / T
                + 1.954263612 * 10 ** 1
                - 2.737830188 * 10 ** -2 * T
                + 1.6261698 * 10 ** -5 * T ** 2
                + 7.0229056 * 10 ** -10 * T ** 3
                - 1.8680009 * 10 ** -13 * T ** 4
                + 2.7150305 * math.log(T)
            )
            return math.exp(Logew) / 100

        elif (
            method == "Buck1981"
        ):  # Buck, A. L., New equations for computing vapor pressure and enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981
            tc = T - 273.15
            return 6.1121 * math.exp(17.502 * tc / (240.97 + tc))

        elif (
            method == "Buck1996"
        ):  # Buck Research Manual (1996); updated equation from Buck1981
            tc = T - 273.15
            return 6.1121 * math.exp((18.678 - tc / 234.5) * tc / (257.14 + tc))

        elif method == "Sonntag":  # (Sonntag, 1994)
            Logew = (
                -6096.9385 / T
                + 16.635794
                - 2.711193 * 10 ** -2 * T
                + 1.673952 * 10 ** -5 * T ** 2
                + 2.433502 * math.log(T)
            )
            return math.exp(Logew)

        elif method == "Magnus-Tetens":  # Murray, 1967
            return 6.1078 * math.exp(17.269388 * (T - 273.16) / (T - 35.86))

        # Bolton, D., The computation of equivalent potential temperature, Monthly Weather Review, 108, 1046-1053, 1980..
        elif method == "Bolton":
            tc = T - 273.15
            return 6.112 * math.exp(17.67 * tc / (tc + 243.5))

        elif method == "Murphy-Koop":  # Murphy and Koop, 2005
            Logew = (
                54.842763
                - 6763.22 / T
                - 4.21 * math.log(T)
                + 0.000367 * T
                + math.tanh(0.0415 * (T - 218.8))
                * (53.878 - 1331.22 / T - 9.44523 * math.log(T) + 0.014025 * T)
            )
            return math.exp(Logew) / 100

        elif (
            method == "IAPWS"
        ):  # IAPWS (International Association for the Properties of Water and Steam) Formulation 1995 (Wagner and Pruß, 2002)
            nu = 1 - T / 647.096
            Log_ew_22064e6 = (
                647.096
                / T
                * (
                    (
                        -7.85951783 * nu
                        + 1.84408259 * nu ** 1.5
                        - 11.7866497 * nu ** 3
                        + 22.6807411 * nu ** 3.5
                        - 15.9618719 * nu ** 4
                        + 1.80122502 * nu ** 7.5
                    )
                )
            )
            return math.exp(Log_ew_22064e6) * 22.064e6 / 100

        else:
            raise ValueError(
                "Method not found for water surface. Support methods: ['Goff-Gratch', 'CIMO', 'WMO2012', 'Hyland-Wexler', 'Buck1981', 'Buck1996', 'Sonntag', 'Magnus-Tetens', 'Bolton', 'Murphy-Koop', 'IAPWS']"
            )

    elif surface == "ice":
        if (
            method == "Goff-Gratch"
        ):  # (Smithsonian Tables, 1984, after Goff and Gratch, 1946)
            Log10ei = (
                -9.09718 * (273.16 / T - 1)
                - 3.56654 * math.log10(273.16 / T)
                + 0.876793 * (1 - T / 273.16)
                + math.log10(6.1071)
            )
            return 10 ** Log10ei

        elif (
            method == "CIMO"
        ):  # Guide to Meteorological Instruments and Methods of Observation (CIMO Guide)
            tc = T - 273.15
            return 6.112 * math.exp(22.46 * tc / (272.62 + tc))

        # Hyland, R. W. and A. Wexler, Formulations for the Thermodynamic Properties of the saturated Phases of H2O from 173.15K to 473.15K, ASHRAE Trans, 89(2A), 500-519, 1983.
        elif method == "Hyland-Wexler":
            Logei = (
                -0.56745359 * 10 ** 4 / T
                + 0.63925247 * 10 ** 1
                - 0.96778430 * 10 ** -2 * T
                + 0.62215701 * 10 ** -6 * T ** 2
                + 0.20747825 * 10 ** -8 * T ** 3
                - 0.94840240 * 10 ** -12 * T ** 4
                + 0.41635019 * 10 ** 1 * math.log(T)
            )
            return math.exp(Logei) / 100

        elif (
            method == "Buck1981"
        ):  # Buck, A. L., New equations for computing vapor pressure and enhancement factor, J. Appl. Meteorol., 20, 1527-1532, 1981
            tc = T - 273.15
            return 6.1115 * math.exp(22.452 * tc / (272.55 + tc))

        elif (
            method == "Buck1996"
        ):  # Buck Research Manual (1996); updated equation from Buck1981
            tc = T - 273.15
            return 6.1115 * math.exp((23.036 - tc / 333.7) * tc / (279.82 + tc))

        elif method == "Sonntag":  # (Sonntag, 1994)
            Logei = (
                -6024.5282 / T
                + 24.721994
                + 1.0613868 * 10 ** -2 * T
                - 1.3198825 * 10 ** -5 * T ** 2
                - 0.49382577 * math.log(T)
            )
            return math.exp(Logei)

        elif method == "Magnus-Tetens":  # Murray, 1967
            return 6.1078 * math.exp(21.8745584 * (T - 273.16) / (T - 7.66))

        elif method == "Murphy-Koop":  # (Murphy and Koop, 2005)
            Logei = 9.550426 - 5723.265 / T + \
                3.53068 * math.log(T) - 0.00728332 * T
            return math.exp(Logei) / 100

        elif method == "Marti-Mauersberger":  # (Marti and Mauersberger, 1993)
            Log10ei = -2663.5 / T + 12.537
            return 10 ** Log10ei / 100

        else:
            raise ValueError(
                "Method not found for ice surface. Support methods: ['Goff-Gratch', 'CIMO', 'Hyland-Wexler', 'Buck1981', 'Buck1996', 'Sonntag', 'Magnus-Tetens', 'Murphy-Koop', 'Marti-Mauersberger']"
            )

    else:
        raise ValueError("surface must be either 'water' or 'ice'")
