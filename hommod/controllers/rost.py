from math import pow, exp


def get_min_identity(nalign):
    """ 
    http://dx.doi.org/10.1093/protein/12.2.85
    Rost curve, same as in HOPE.
    """
    if nalign <= 0:
        return float('inf')

    n = float(nalign)
    return 480 * pow(n, -0.32 * (1 + exp(-n / 1000)))
