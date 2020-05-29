import numpy as np
from mvg_moments import *

def collision_rv_moments_six(input_moments, input_deterministic):
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
    xPow10 = input_moments[(10, 0)]
    xPow9_yPow1 = input_moments[(9, 1)]
    xPow8_yPow2 = input_moments[(8, 2)]
    xPow7_yPow3 = input_moments[(7, 3)]
    xPow6_yPow4 = input_moments[(6, 4)]
    xPow5_yPow5 = input_moments[(5, 5)]
    xPow4_yPow6 = input_moments[(4, 6)]
    xPow3_yPow7 = input_moments[(3, 7)]
    xPow2_yPow8 = input_moments[(2, 8)]
    xPow1_yPow9 = input_moments[(1, 9)]
    yPow10 = input_moments[(0, 10)]
    xPow12 = input_moments[(12, 0)]
    xPow11_yPow1 = input_moments[(11, 1)]
    xPow10_yPow2 = input_moments[(10, 2)]
    xPow9_yPow3 = input_moments[(9, 3)]
    xPow8_yPow4 = input_moments[(8, 4)]
    xPow7_yPow5 = input_moments[(7, 5)]
    xPow6_yPow6 = input_moments[(6, 6)]
    xPow5_yPow7 = input_moments[(5, 7)]
    xPow4_yPow8 = input_moments[(4, 8)]
    xPow3_yPow9 = input_moments[(3, 9)]
    xPow2_yPow10 = input_moments[(2, 10)]
    xPow1_yPow11 = input_moments[(1, 11)]
    yPow12 = input_moments[(0, 12)]
    q11 = input_deterministic["q11"]
    q12 = input_deterministic["q12"]
    q22 = input_deterministic["q22"]

    moments = dict()

    # Moment expressions.
    moments[1] = q11*xPow2 + 2*q12*xPow1_yPow1 + q22*yPow2 - 1.0
    moments[2] = q11**2*xPow4 + 4*q11*q12*xPow3_yPow1 - 2*q11*xPow2 + 4*q12*q22*xPow1_yPow3 - 4*q12*xPow1_yPow1 + q22**2*yPow4 - 2*q22*yPow2 + xPow2_yPow2*(2*q11*q22 + 4*q12**2) + 1.0
    moments[3] = q11**3*xPow6 + 6*q11**2*q12*xPow5_yPow1 - 3*q11**2*xPow4 - 12*q11*q12*xPow3_yPow1 + 3*q11*xPow2 + 6*q12*q22**2*xPow1_yPow5 - 12*q12*q22*xPow1_yPow3 + 6*q12*xPow1_yPow1 + q22**3*yPow6 - 3*q22**2*yPow4 + 3*q22*yPow2 + xPow2_yPow2*(-6*q11*q22 - 12*q12**2) + xPow2_yPow4*(3*q11*q22**2 + 12*q12**2*q22) + xPow3_yPow3*(12*q11*q12*q22 + 8*q12**3) + xPow4_yPow2*(3*q11**2*q22 + 12*q11*q12**2) - 1.0
    moments[4] = q11**4*xPow8 + 8*q11**3*q12*xPow7_yPow1 - 4*q11**3*xPow6 - 24*q11**2*q12*xPow5_yPow1 + 6*q11**2*xPow4 + 24*q11*q12*xPow3_yPow1 - 4*q11*xPow2 + 8*q12*q22**3*xPow1_yPow7 - 24*q12*q22**2*xPow1_yPow5 + 24*q12*q22*xPow1_yPow3 - 8*q12*xPow1_yPow1 + q22**4*yPow8 - 4*q22**3*yPow6 + 6*q22**2*yPow4 - 4*q22*yPow2 + xPow2_yPow2*(12*q11*q22 + 24*q12**2) + xPow2_yPow4*(-12*q11*q22**2 - 48*q12**2*q22) + xPow2_yPow6*(4*q11*q22**3 + 24*q12**2*q22**2) + xPow3_yPow3*(-48*q11*q12*q22 - 32*q12**3) + xPow3_yPow5*(24*q11*q12*q22**2 + 32*q12**3*q22) + xPow4_yPow2*(-12*q11**2*q22 - 48*q11*q12**2) + xPow4_yPow4*(6*q11**2*q22**2 + 48*q11*q12**2*q22 + 16*q12**4) + xPow5_yPow3*(24*q11**2*q12*q22 + 32*q11*q12**3) + xPow6_yPow2*(4*q11**3*q22 + 24*q11**2*q12**2) + 1.0
    moments[5] = q11**5*xPow10 + 10*q11**4*q12*xPow9_yPow1 - 5*q11**4*xPow8 - 40*q11**3*q12*xPow7_yPow1 + 10*q11**3*xPow6 + 60*q11**2*q12*xPow5_yPow1 - 10*q11**2*xPow4 - 40*q11*q12*xPow3_yPow1 + 5*q11*xPow2 + 10*q12*q22**4*xPow1_yPow9 - 40*q12*q22**3*xPow1_yPow7 + 60*q12*q22**2*xPow1_yPow5 - 40*q12*q22*xPow1_yPow3 + 10*q12*xPow1_yPow1 + q22**5*yPow10 - 5*q22**4*yPow8 + 10*q22**3*yPow6 - 10*q22**2*yPow4 + 5*q22*yPow2 + xPow2_yPow2*(-20*q11*q22 - 40*q12**2) + xPow2_yPow4*(30*q11*q22**2 + 120*q12**2*q22) + xPow2_yPow6*(-20*q11*q22**3 - 120*q12**2*q22**2) + xPow2_yPow8*(5*q11*q22**4 + 40*q12**2*q22**3) + xPow3_yPow3*(120*q11*q12*q22 + 80*q12**3) + xPow3_yPow5*(-120*q11*q12*q22**2 - 160*q12**3*q22) + xPow3_yPow7*(40*q11*q12*q22**3 + 80*q12**3*q22**2) + xPow4_yPow2*(30*q11**2*q22 + 120*q11*q12**2) + xPow4_yPow4*(-30*q11**2*q22**2 - 240*q11*q12**2*q22 - 80*q12**4) + xPow4_yPow6*(10*q11**2*q22**3 + 120*q11*q12**2*q22**2 + 80*q12**4*q22) + xPow5_yPow3*(-120*q11**2*q12*q22 - 160*q11*q12**3) + xPow5_yPow5*(60*q11**2*q12*q22**2 + 160*q11*q12**3*q22 + 32*q12**5) + xPow6_yPow2*(-20*q11**3*q22 - 120*q11**2*q12**2) + xPow6_yPow4*(10*q11**3*q22**2 + 120*q11**2*q12**2*q22 + 80*q11*q12**4) + xPow7_yPow3*(40*q11**3*q12*q22 + 80*q11**2*q12**3) + xPow8_yPow2*(5*q11**4*q22 + 40*q11**3*q12**2) - 1.0
    moments[6] = q11**6*xPow12 + 12*q11**5*q12*xPow11_yPow1 - 6*q11**5*xPow10 - 60*q11**4*q12*xPow9_yPow1 + 15*q11**4*xPow8 + 120*q11**3*q12*xPow7_yPow1 - 20*q11**3*xPow6 - 120*q11**2*q12*xPow5_yPow1 + 15*q11**2*xPow4 + 60*q11*q12*xPow3_yPow1 - 6*q11*xPow2 + 12*q12*q22**5*xPow1_yPow11 - 60*q12*q22**4*xPow1_yPow9 + 120*q12*q22**3*xPow1_yPow7 - 120*q12*q22**2*xPow1_yPow5 + 60*q12*q22*xPow1_yPow3 - 12*q12*xPow1_yPow1 + q22**6*yPow12 - 6*q22**5*yPow10 + 15*q22**4*yPow8 - 20*q22**3*yPow6 + 15*q22**2*yPow4 - 6*q22*yPow2 + xPow10_yPow2*(6*q11**5*q22 + 60*q11**4*q12**2) + xPow2_yPow10*(6*q11*q22**5 + 60*q12**2*q22**4) + xPow2_yPow2*(30*q11*q22 + 60*q12**2) + xPow2_yPow4*(-60*q11*q22**2 - 240*q12**2*q22) + xPow2_yPow6*(60*q11*q22**3 + 360*q12**2*q22**2) + xPow2_yPow8*(-30*q11*q22**4 - 240*q12**2*q22**3) + xPow3_yPow3*(-240*q11*q12*q22 - 160*q12**3) + xPow3_yPow5*(360*q11*q12*q22**2 + 480*q12**3*q22) + xPow3_yPow7*(-240*q11*q12*q22**3 - 480*q12**3*q22**2) + xPow3_yPow9*(60*q11*q12*q22**4 + 160*q12**3*q22**3) + xPow4_yPow2*(-60*q11**2*q22 - 240*q11*q12**2) + xPow4_yPow4*(90*q11**2*q22**2 + 720*q11*q12**2*q22 + 240*q12**4) + xPow4_yPow6*(-60*q11**2*q22**3 - 720*q11*q12**2*q22**2 - 480*q12**4*q22) + xPow4_yPow8*(15*q11**2*q22**4 + 240*q11*q12**2*q22**3 + 240*q12**4*q22**2) + xPow5_yPow3*(360*q11**2*q12*q22 + 480*q11*q12**3) + xPow5_yPow5*(-360*q11**2*q12*q22**2 - 960*q11*q12**3*q22 - 192*q12**5) + xPow5_yPow7*(120*q11**2*q12*q22**3 + 480*q11*q12**3*q22**2 + 192*q12**5*q22) + xPow6_yPow2*(60*q11**3*q22 + 360*q11**2*q12**2) + xPow6_yPow4*(-60*q11**3*q22**2 - 720*q11**2*q12**2*q22 - 480*q11*q12**4) + xPow6_yPow6*(20*q11**3*q22**3 + 360*q11**2*q12**2*q22**2 + 480*q11*q12**4*q22 + 64*q12**6) + xPow7_yPow3*(-240*q11**3*q12*q22 - 480*q11**2*q12**3) + xPow7_yPow5*(120*q11**3*q12*q22**2 + 480*q11**2*q12**3*q22 + 192*q11*q12**5) + xPow8_yPow2*(-30*q11**4*q22 - 240*q11**3*q12**2) + xPow8_yPow4*(15*q11**4*q22**2 + 240*q11**3*q12**2*q22 + 240*q11**2*q12**4) + xPow9_yPow3*(60*q11**4*q12*q22 + 160*q11**3*q12**3) + 1.0
    return moments


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


