import scipy.io as sio
from collision_rv_moments import *
import numpy as np
from utils import *

def generate_data():
    # Specify the Q matrix.
    Q = np.array([[1.0, 0.0], [0.0, 1.0]])
    input_deterministic = dict()
    input_deterministic["q11"] = Q[0, 0]
    input_deterministic["q12"] = Q[0, 1]
    input_deterministic["q22"] = Q[1, 1]


    mu = np.array([2.0, 0.0])
    sigma = np.array([[0.5, 0.0], [0.0, 0.5]])
    ys = np.linspace(-10, 10, 30)
    collision_moments = []
    for y in ys:
        mu[1] = y
        input_moments = mvg_moments(mu, sigma, 12) # need 2 * n momnets of the MVG
        out = collision_rv_moments_six(input_moments, input_deterministic)
        collision_m = offset_moments(out, -1)
        collision_moments.append(list(collision_m.values()))
    print(collision_moments)
    sio.savemat('test.mat', {'collision_moments': collision_moments})

generate_data()