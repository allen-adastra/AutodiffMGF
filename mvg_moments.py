import casadi
import numpy as np

def constant_sum_tuples(length, total_sum):
    """ Generate all tuples that sum to a given value.

    Args:
        length ([type]): [description]
        total_sum ([type]): [description]

    Yields:
        [type]: [description]
    """
    if length == 1:
        yield (total_sum,)
    else:
        for value in range(total_sum + 1):
            for permutation in constant_sum_tuples(length - 1, total_sum - value):
                yield (value,) + permutation

def mvg_mgf(t, mu, sigma):
    """Multivariate Gaussian Moment Generating Function.

    Args:
        t (1D array): 
        mu (1D array): 
        sigma (2D array): 
    """
    return casadi.exp(mu.T @ t + 0.5 * t.T @ sigma @ t)

def tensor_to_multi_idx(tensor_idx, dimension):
    return tuple(tensor_idx.count(i) for i in range(dimension))

def multi_to_tensor_idx(multi_idx, dimension):
    # Build the index of the corresponding element in the moment tensor.
    tensor_idx = tuple()
    for var_number in range(dimension):
        # var_number corresponds to an element of the random vector.
        # The order of the moment of "var" corresponds to the number of times 
        # the index of "var" shows up in the index of the moment matrix.
        tensor_idx += multi_idx[var_number]*(var_number,)
    return tensor_idx

def mvg_moment_array_functions(max_order, dimension):
    """ 
    Args:
        max_order ([type]): [description]
        dimension ([type]): [description]

    Returns:
        [type]: Resulting functions have signature (t, mu, sigma) where 
                t and mu are vectors, and sigma is the covariance matrix.
    """
    # Declare variables.
    mu = casadi.MX.sym("mu", dimension, 1)
    sigma = casadi.MX.sym("sigma", dimension, dimension)
    t = casadi.MX.sym('t', dimension, 1)

    # This dictionary stores our resulting functions.
    moment_array_functions = dict()

    # Start auto-differentiating the multivariate Gaussian MGF>
    foo = mvg_mgf(t, mu, sigma)
    for i in range(max_order):
        order = i+1
        foo = casadi.jacobian(foo, t) #TODO: current matrix output form doesn't quite make sense
        fun = casadi.Function('moment_array_order' + str(order), [t, mu, sigma], [foo])
        moment_array_functions[order] = fun
    return moment_array_functions

def mvg_moment_tensors(mu, sigma, max_order, dimension):
    """

    Args:
        mu ([type]): [description]
        sigma ([type]): [description]
        max_order ([type]): [description]
        dimension ([type]): [description]

    Returns:
        [type]: [description]
    """
    moment_array_funcs = mvg_moment_array_functions(max_order, dimension)
    t = np.zeros(dimension)
    moment_tensors = dict()
    for order, fun in moment_array_funcs.items():
        moment_tensors[order] = np.array(fun(t, mu, sigma)).reshape(order*(2,))
    return moment_tensors

def mvg_moments(mu, sigma, max_order):
    """ Compute all moments up to the specified order. Returns moments
        as a dictionary with moment multi-index as the key.

    Args:
        mu ([type]): [description]
        sigma ([type]): [description]
        max_order ([type]): [description]
    """
    # Compute the moment tensors.
    dimension = mu.size
    if sigma.shape != (dimension, dimension):
        raise Exception("Incompatible dimensions for mu and sigma.")
    moment_tensors = mvg_moment_tensors(mu, sigma, max_order, dimension)

    # Populate moment_values with moments up to "max_order".
    moment_values = dict()
    for order in range(1, max_order+1):
        # Retrieve all moments of order "order" from the moment tensor.
        moment_tensor = moment_tensors[order]

        # All possible moment multi-indices of the current order.
        possible_multi_idxs = constant_sum_tuples(dimension, order)
        
        for multi_idx in possible_multi_idxs:
            # Retrieve the desired moment value by indexing into the moment tensor.
            tensor_idx = multi_to_tensor_idx(multi_idx, dimension)
            moment_values[multi_idx] = moment_tensor[tensor_idx]
    return moment_values