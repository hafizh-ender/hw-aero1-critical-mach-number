from numpy import sqrt


####################################################################################################
#                                                                                                  #
#                                  Functions and Constants                                         #
#                                                                                                  #
####################################################################################################

defaultCp0 = -1.0
defaultGamma = 1.4

minimum_float = 1.175494351e-38             # used to avoid division by zero
toleratedDifferenceFromZero = 10e-16                

MmaxPrandtlGlauert = 1                      # location of Prandtl-Glauert rule asymptote

def MmaxLaitone(Cp0=-defaultCp0, gamma=defaultGamma):
    # Get the location of Laitone rule asymptote
    return sqrt((2-Cp0-sqrt((2-Cp0)**2-4*((gamma-1)*Cp0)))/((gamma-1)*Cp0))

def MmaxKarmanTsien(Cp0=-defaultCp0):
    # Get the location of Karman-Tsien rule asymptote
    return sqrt(1 - ((1-sqrt(1+Cp0*(Cp0-2)))/(Cp0-2))**2)

def correction(M):
    return sqrt(1 - M**2)

def criticalPressureCoefficients(Mcr, gamma_=defaultGamma):
    return 2 / (gamma_ * Mcr**2) * (((1 + (gamma_ - 1) / 2 * Mcr**2) / (1 + (gamma_ - 1) / 2))**(gamma_ / (gamma_ - 1)) - 1)

def prandtlGlauert(Minf, Cp0_=defaultCp0):
    return Cp0_ / correction(Minf)

def laitone(Minf, Cp0_=defaultCp0, gamma_=defaultGamma):
    return Cp0_ / (correction(Minf) + (Minf**2 * (1 + (gamma_ - 1)/2 * Minf**2) / (2 * correction(Minf))) * Cp0_)

def karmanTsien(Minf, Cp0_=defaultCp0):
    return Cp0_ / (correction(Minf) + (Minf**2 / (1 + correction(Minf))) * Cp0_ / 2)

def rootFuncPrandtlGlauert(Mcr, Cp0_=defaultCp0, gamma_=defaultGamma):
    return criticalPressureCoefficients(Mcr, gamma_) - prandtlGlauert(Mcr, Cp0_)

def rootFuncLaitone(Mcr, Cp0_=defaultCp0, gamma_=defaultGamma):
    return criticalPressureCoefficients(Mcr, gamma_) - laitone(Mcr, Cp0_, gamma_)

def rootFuncKarmanTsien(Mcr, Cp0_=defaultCp0, gamma_=defaultGamma):
    return criticalPressureCoefficients(Mcr, gamma_) - karmanTsien(Mcr, Cp0_)

def bisectionMethodPrandtlGlauert(Mil, Miu, iter, Cp0_=defaultCp0, gamma_=defaultGamma):
    # Get the initial guess
    Mil_ = Mil
    Miu_ = Miu
    Mi_ = (Mil_ + Miu_) / 2

    # Get the number of iterations
    n_ = iter
    i = n_

    # Iterate until maximum number of iterations is reached or the tolerance is met
    while (i > 0 and abs(rootFuncPrandtlGlauert(Mi_, Cp0_, gamma_)) > toleratedDifferenceFromZero):
        if (rootFuncPrandtlGlauert(Mil_, Cp0_, gamma_) * rootFuncPrandtlGlauert(Mi_, Cp0_, gamma_) < 0):
            Miu_ = Mi_
        else:
            Mil_ = Mi_
        
        Mi_ = (Mil_ + Miu_) / 2

        i = i - 1

    return Mi_

def bisectionMethodLaitone(Mil, Miu, iter, Cp0_=defaultCp0, gamma_=defaultGamma):
    # Get the initial guess
    Mil_ = Mil
    Miu_ = Miu
    Mi_ = (Mil_ + Miu_) / 2

    # Get the number of iterations
    n_ = iter
    i = n_

    # Iterate until maximum number of iterations is reached or the tolerance is met
    while (i > 0 and abs(rootFuncLaitone(Mi_, Cp0_, gamma_)) > toleratedDifferenceFromZero):
        if (rootFuncLaitone(Mil_, Cp0_, gamma_) * rootFuncLaitone(Mi_, Cp0_, gamma_) < 0):
            Miu_ = Mi_
        else:
            Mil_ = Mi_
        
        Mi_ = (Mil_ + Miu_) / 2

        i = i - 1

    return Mi_

def bisectionMethodKarmanTsien(Mil, Miu, iter, Cp0_=defaultCp0, gamma_=defaultGamma):
    # Get the initial guess
    Mil_ = Mil
    Miu_ = Miu
    Mi_ = (Mil_ + Miu_) / 2

    # Get the number of iterations
    n_ = iter
    i = n_

    # Iterate until maximum number of iterations is reached or the tolerance is met
    while (i > 0 and abs(rootFuncKarmanTsien(Mi_, Cp0_, gamma_)) > toleratedDifferenceFromZero):
        if (rootFuncKarmanTsien(Mil_, Cp0_, gamma_) * rootFuncKarmanTsien(Mi_, Cp0_, gamma_) < 0):
            Miu_ = Mi_
        else:
            Mil_ = Mi_
        
        Mi_ = (Mil_ + Miu_) / 2

        i = i - 1

    return Mi_
