from mvg_moments import *
from collision_rv_moments import *

from utils import *


def test_moment_tensor_consistent(tolerance=1e-10):
    """ Test that elements of the moment tensor are consistent. That is,
        for any two moment tensor indices, if they have the same corresponding
        multi-index, then they have the same values.
    """
    # Test parameters.
    max_order = 8
    dimension = 2
    mu = np.array([-2.0, 1.0])
    sigma = np.array([[1.0, 0.1], [0.1, 0.5]])

    moment_tensors = mvg_moment_tensors(mu, sigma, max_order, dimension)
    moment_values = mvg_moments(mu, sigma, max_order)
    for moment_tensor in moment_tensors.values():
        for tensor_idx in np.ndindex(moment_tensor.shape):
            multi_idx = tensor_to_multi_idx(tensor_idx, dimension)
            assert abs(moment_tensor[tensor_idx]-moment_values[multi_idx])<=tolerance

def test_collision_rv_moments(tolerance = 1e-10):
    """ Test the collision rv moments that are computed against
        a special case where the distribution is chi-squared.
    """
    # Specify the MVG x and Q s.t. Q(x) is a 2dof chi-squared 
    mu = np.array([0.0, 0.0])
    sigma = np.array([[1.0, 0.0], [0.0, 1.0]])
    Q = np.array([[1.0, 0.0], [0.0, 1.0]])
    input_deterministic = dict()
    input_deterministic["q11"] = Q[0, 0]
    input_deterministic["q12"] = Q[0, 1]
    input_deterministic["q22"] = Q[1, 1]

    # Test the function to compute moments of Q(x) - 1.
    input_moments = mvg_moments(mu, sigma, 8)
    computed_moments = collision_rv_moments(input_moments, input_deterministic)

    # Compute ground truth.
    Qmoments = {order : chi_square_moments(order, 2) for order in range(1, 5)}
    ground_truth = offset_moments(Qmoments, -1)

    # Check results.
    for i in range(1, 5):
        assert abs(ground_truth[i] - computed_moments[i]) <= tolerance