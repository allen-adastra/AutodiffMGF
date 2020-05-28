import numpy as np
from mvg_moments import *

def collision_rv_moments(input_moments, input_deterministic):
    # Parse required inputs.
    xPow2 = input_moments[(2, 0)]
    xPow1_yPow1 = input_moments[(1, 1)]
    yPow2 = input_moments[(0, 2)]
    xPow4 = input_moments[(4, 0)]
    xPow3_yPow1 = input_moments[(3, 1)]
    xPow2_yPow2 = input_moments[(2, 2)]
    xPow1_yPow3 = input_moments[(1, 3)]
    yPow4 = input_moments[(0, 4)]
    xPow6 = input_moments[(6, 0)]
    xPow5_yPow1 = input_moments[(5, 1)]
    xPow4_yPow2 = input_moments[(4, 2)]
    xPow3_yPow3 = input_moments[(3, 3)]
    xPow2_yPow4 = input_moments[(2, 4)]
    xPow1_yPow5 = input_moments[(1, 5)]
    yPow6 = input_moments[(0, 6)]
    xPow8 = input_moments[(8, 0)]
    xPow7_yPow1 = input_moments[(7, 1)]
    xPow6_yPow2 = input_moments[(6, 2)]
    xPow5_yPow3 = input_moments[(5, 3)]
    xPow4_yPow4 = input_moments[(4, 4)]
    xPow3_yPow5 = input_moments[(3, 5)]
    xPow2_yPow6 = input_moments[(2, 6)]
    xPow1_yPow7 = input_moments[(1, 7)]
    yPow8 = input_moments[(0, 8)]

    q11 = input_deterministic["q11"]
    q12 = input_deterministic["q12"]
    q22 = input_deterministic["q22"]
    
    moments = dict()

    # Moment expressions.
    moments[1] = q11*xPow2 + 2*q12*xPow1_yPow1 + q22*yPow2 - 1.0
    moments[2] = q11**2*xPow4 + 4*q11*q12*xPow3_yPow1 - 2*q11*xPow2 + 4*q12*q22*xPow1_yPow3 - 4*q12*xPow1_yPow1 + q22**2*yPow4 - 2*q22*yPow2 + xPow2_yPow2*(2*q11*q22 + 4*q12**2) + 1.0
    moments[3] = q11**3*xPow6 + 6*q11**2*q12*xPow5_yPow1 - 3*q11**2*xPow4 - 12*q11*q12*xPow3_yPow1 + 3*q11*xPow2 + 6*q12*q22**2*xPow1_yPow5 - 12*q12*q22*xPow1_yPow3 + 6*q12*xPow1_yPow1 + q22**3*yPow6 - 3*q22**2*yPow4 + 3*q22*yPow2 + xPow2_yPow2*(-6*q11*q22 - 12*q12**2) + xPow2_yPow4*(3*q11*q22**2 + 12*q12**2*q22) + xPow3_yPow3*(12*q11*q12*q22 + 8*q12**3) + xPow4_yPow2*(3*q11**2*q22 + 12*q11*q12**2) - 1.0
    moments[4] = q11**4*xPow8 + 8*q11**3*q12*xPow7_yPow1 - 4*q11**3*xPow6 - 24*q11**2*q12*xPow5_yPow1 + 6*q11**2*xPow4 + 24*q11*q12*xPow3_yPow1 - 4*q11*xPow2 + 8*q12*q22**3*xPow1_yPow7 - 24*q12*q22**2*xPow1_yPow5 + 24*q12*q22*xPow1_yPow3 - 8*q12*xPow1_yPow1 + q22**4*yPow8 - 4*q22**3*yPow6 + 6*q22**2*yPow4 - 4*q22*yPow2 + xPow2_yPow2*(12*q11*q22 + 24*q12**2) + xPow2_yPow4*(-12*q11*q22**2 - 48*q12**2*q22) + xPow2_yPow6*(4*q11*q22**3 + 24*q12**2*q22**2) + xPow3_yPow3*(-48*q11*q12*q22 - 32*q12**3) + xPow3_yPow5*(24*q11*q12*q22**2 + 32*q12**3*q22) + xPow4_yPow2*(-12*q11**2*q22 - 48*q11*q12**2) + xPow4_yPow4*(6*q11**2*q22**2 + 48*q11*q12**2*q22 + 16*q12**4) + xPow5_yPow3*(24*q11**2*q12*q22 + 32*q11*q12**3) + xPow6_yPow2*(4*q11**3*q22 + 24*q11**2*q12**2) + 1.0
    return moments


