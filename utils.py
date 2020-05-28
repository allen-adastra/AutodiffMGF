import math
import sympy as sp

def chi_square_moments(order, dof):
    # https://en.wikipedia.org/wiki/Chi-square_distribution#Noncentral_moments
    num = 2**order * math.gamma(order + 0.5*dof)
    denom = math.gamma(0.5*dof)
    return num/denom

def offset_moments(moments, offset):
    max_order = max(moments.keys())
    m = sp.Symbol("m")
    offset_moments = dict()
    for order in range(1, max_order+1):
        offset_expression = sp.poly((m + offset)**order, m)
        # Substitute moments accordingly.
        for i in range(order, 0, -1):
            offset_expression = offset_expression.subs(m**i, moments[i])
        offset_moments[order] = offset_expression
    return offset_moments