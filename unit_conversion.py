import numpy as np

# Constants
R = 8.3144621 # J/(mol*K)
T = 298.15 # K

def delta_G_to_pKd(delta_G):
    """ delta_G in kcal/mol, pKd in M
    """
    return -delta_G * 1000 / (R * T * np.log(10)) * 4.184

def pKd_to_delta_G(pKd):
    """ pKd in M, delta_G in kcal/mol
    """
    return -pKd * R * T * np.log(10) / 1000 / 4.184

if __name__ == '__main__':
    # print(pKd_to_delta_G(9))
    # print(delta_G_to_pKd(pKd_to_delta_G(9)))

    print(pKd_to_delta_G(0.04))