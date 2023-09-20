"""
Greeks of Binomial Tree Models
"""


class greekType():
    def Delta(S, r, vol, trade, strike, n, calib):
        delta = (binomialPricer(S+0.001*S, r, vol, trade, n, calib) 
                - binomialPricer(S-0.001*S, r, vol, trade, n, calib)
                ) / (2*0.001*S)
        return delta

    def Gamma(S, r, vol, trade, strike, n, calib):
        gamma = (binomialPricer(S+0.001*S, r, vol, trade, n, calib) 
                - 2 * binomialPricer(S, r, vol, trade, n, calib) 
                + binomialPricer(S-0.001*S, r, vol, trade, n, calib)
                ) / (0.001*S)**2
        return gamma

    def Vega(S, r, vol, trade, strike, n, calib):
        vega = (binomialPricer(S, r, vol+0.001, trade, n, calib)
                - binomialPricer(S, r, vol-0.001, trade, n, calib)
                ) / (0.001*2)
        return vega

    def Theta(S, r, vol, trade, strike, n, calib):
        Theta = (binomialPricer1(S, r, vol, trade, n, calib)
                - binomialPricer(S, r, vol, trade, n, calib)
                ) / (0.004*2)
        return Theta

    def Rho(S, r, vol, trade, strike, n, calib):
        rho = (binomialPricer(S, r+0.001, vol, trade, n, calib)
                - binomialPricer(S, r-0.001, vol, trade, n, calib)
                    )/(2*0.001)
        return rho


def binomialGreeks(S, r, vol, T, strike, greekType) -> float:
    ks = range(50,150,1)
    greeks_list = []
    
    for trade in [AmericanOption, EuropeanOption]:
        for payoffType in [PayoffType.Call, PayoffType.Put]:
            for calib in [crrCalib, jrrnCalib, jreqCalib, tianCalib]:
                greek = [greekType(S, r, vol, trade(T, k, payoffType), k, n, calib) for k in ks]
                greeks_list.append(greek)
    
    return greeks_list