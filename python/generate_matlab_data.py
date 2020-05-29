import scipy.io as sio
from collision_rv_moments import *
import numpy as np
from utils import *
import math
import matplotlib.pyplot as plt
from risk_assess.random_objects.random_variables import MultivariateNormal
from risk_assess.random_objects.quad_forms import MvnQuadForm

def rotation_matrix(theta):
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[c, -s], [s, c]])

def generate_data():
    # Specify the Q matrix.
    Q = np.array([[0.25, 0.0], [0.0, 0.5]])
    input_deterministic = dict()
    input_deterministic["q11"] = Q[0, 0]
    input_deterministic["q12"] = Q[0, 1]
    input_deterministic["q22"] = Q[1, 1]

    dt = 0.1
    x = 3.0
    y = -1.0
    v = 7.0
    theta = 0.3 * math.pi
    sigma = np.array([[0.1, 0.02], [0.02, 0.1]])
    sigma_dot = 0.1
    cross_sigma_dot = 0.02
    accel = 0.0
    steer = 1.25
    steps = 30

    ground_truths = []
    collision_moments = []
    xs = []
    ys = []

    # Ego speed
    ego_drift = 3.3
    for i in range(steps):
        # Update state
        x += dt * v * math.cos(theta)
        y += dt * v * math.sin(theta) - dt * ego_drift
        xs.append(x)
        ys.append(y)
        v += dt * accel
        theta += dt * steer

        # Update mu and sigma
        mu = np.array([x, y])
        sigma += dt * sigma_dot * np.eye(2)
        sigma += dt * cross_sigma_dot * (np.ones((2, 2)) - np.eye(2))
        R = rotation_matrix(-dt * steer)
        sigma = R @ sigma @ R.T

        # Compute ground truth.
        mvn = MultivariateNormal(mu, sigma)
        mvnqf = MvnQuadForm(Q, mvn)
        ground_truths.append(1.0 - mvnqf.upper_tail_probability_imhof(1.0, 1e-8, 1e-8))

        # Compute moments.
        input_moments = mvg_moments(mu, sigma, 12) # need 2 * n momnets of the MVG
        out = collision_rv_moments_six(input_moments, input_deterministic)
        collision_m = offset_moments(out, -1)
        collision_moments.append(list(collision_m.values()))
    sio.savemat('test.mat', {'collision_moments': collision_moments, 'ground_truths': ground_truths})
generate_data()