from ifa_test import ifa_simulation
import logging
import numpy as np
from scipy.optimize import minimize

# Function to evaluate the antenna performance
def evaluation_fun(x, variable_names, fixed_params):
    # Update the variables to be optimized
    params = fixed_params.copy()
    for i, var_name in enumerate(variable_names):
        params[var_name] = x[i]

    # Extract parameters
    ifa_h = params['ifa_h']
    ifa_l = params['ifa_l']
    ifa_fp = params['ifa_fp']
    ifa_w = params['ifa_w']

    # Constants for the simulation
    Sim_CSX = 'IFA.xml'
    showCad = False
    post_proc_only = False
    unit = 1e-3
    substrate_width = 21
    substrate_length = 83.15
    substrate_thickness = 1.5
    gndplane_position = 0
    substrate_cells = 4
    ifa_e = 0.5
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 2.4e9
    center_freq = 2.45e9
    max_freq = 2.5e9
    min_size = 0.2  # Minimum automesh size
    plot = False

    try:
        freq, s11_dB, Zin, P_in, hash_prefix = ifa_simulation(
            Sim_CSX=Sim_CSX,
            showCad=showCad,
            post_proc_only=post_proc_only,
            unit=unit,
            substrate_width=substrate_width,
            substrate_length=substrate_length,
            substrate_thickness=substrate_thickness,
            gndplane_position=gndplane_position,
            substrate_cells=substrate_cells,
            ifa_h=ifa_h,
            ifa_l=ifa_l,
            ifa_fp=ifa_fp,
            ifa_w1=ifa_w,
            ifa_w2=ifa_w,
            ifa_wf=ifa_w,
            ifa_e=ifa_e,
            substrate_epsR=substrate_epsR,
            feed_R=feed_R,
            min_freq=min_freq,
            center_freq=center_freq,
            max_freq=max_freq,
            plot=plot,
            min_size=min_size
        )

        # Get the S11 value at the center frequency
        idx = np.argmin(np.abs(freq - center_freq))
        s11_value = s11_dB[idx]
        # Round value to 2 decimal places
        s11_value = round(s11_value, 2)
        # Log parameters and objective function values
        logging.info(f"ifa_l: {ifa_l}, ifa_h: {ifa_h}, ifa_fp: {ifa_fp}, ifa_w: {ifa_w}, S11 at center frequency: {s11_value} - id {hash_prefix}")
        # Return the magnitude (since we want to minimize it)
        return s11_value

    except Exception as e:
        logging.warning(f"Exception during simulation: {e}.")
        return np.inf  # Return a high value to indicate failure

def evaluation_fun_scalar(x_scalar, variable_names, fixed_params):
    return evaluation_fun([x_scalar], variable_names, fixed_params)

if __name__ == "__main__":
    # Setup file logging
    logging.basicConfig(filename='halfsearch_optimization_log.txt', level=logging.INFO, format='%(asctime)s - %(message)s')

    # Fixed parameters
    fixed_params = {
        'ifa_l': 14.0,  # Initial value
        'ifa_h': 7.5,  
        'ifa_fp': 5.5,
        'ifa_w': 0.7
    }

    # Define bounds for each variable you want to optimize
    variable_bounds = {
        'ifa_l': (10.0, 19.5),
        'ifa_h': (1.0, 16.0),
        'ifa_fp': (1.0, 10.0),
        'ifa_w': (0.4, 1.0)
    }

    # Choose variables to optimize
    variable_names = ['ifa_l', 'ifa_h', 'ifa_fp', 'ifa_w']  # List variables you want to optimize
    bcplist = variable_names

    for var_name in bcplist:
        logging.info(f"Optimizing {var_name} within bounds: {variable_bounds[var_name]}")
        variable_names = [var_name]
        bounds = [variable_bounds[var_name]]
        fixed_params_copy = fixed_params.copy()

        # Generate 5 points along the bounds for the variable
        lower_bound, upper_bound = variable_bounds[var_name]
        num_points = 3
        test_points = np.linspace(lower_bound, upper_bound, num_points)

        # Evaluate the objective function at each point
        objective_values = []
        for point in test_points:
            x = [point]
            obj_value = evaluation_fun(x, variable_names, fixed_params_copy)
            objective_values.append(obj_value)
            logging.info(f"Testing {var_name}={point}, Objective function value: {obj_value}")

        # Find the point with the lowest objective function value
        min_index = np.argmin(objective_values)
        best_initial_value = test_points[min_index]
        logging.info(f"Best initial value for {var_name}: {best_initial_value}, Objective function value: {objective_values[min_index]}")

        # Optimization function using half-search method
        # Set initial bounds for optimization
        lower_bound, upper_bound = variable_bounds[var_name]

        # Perform scalar minimization within the bounds
        from scipy.optimize import minimize_scalar

        result = minimize_scalar(
            evaluation_fun_scalar,
            bounds=(lower_bound, upper_bound),
            method='bounded',
            args=(variable_names, fixed_params_copy),
            options={'xatol': 0.1}  # Stop when the interval is smaller than 0.1mm
        )

        fixed_params[var_name] = result.x
        logging.info(f"Optimized {var_name}: {fixed_params[var_name]}, Objective function value: {result.fun}")

    # Log the final optimized parameters
    logging.info("Final optimized parameters:")
    for var_name in bcplist:
        logging.info(f"{var_name}: {fixed_params[var_name]}")
    logging.info(f"Final objective function value: {result.fun}")

    # Print the result to the console
    print("Optimization Result:")
    for var_name in bcplist:
        print(f"{var_name}: {fixed_params[var_name]}")
    print(f"Objective function value (S11 at center frequency): {result.fun}")
