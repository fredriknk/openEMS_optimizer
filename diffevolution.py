from ifa_test import ifa_simulation
import logging
import numpy as np
from scipy.optimize import differential_evolution

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
        freq, s11_dB, Zin, P_in = ifa_simulation(
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
        s11_value = round(s11_value, 4)
        # Log parameters and objective function value
        logging.info(f"ifa_l: {ifa_l}, ifa_h: {ifa_h}, ifa_fp: {ifa_fp}, ifa_w: {ifa_w}, S11 at center frequency: {s11_value}")
        # Return the magnitude (since we want to minimize it)
        return s11_value

    except Exception as e:
        logging.warning(f"Exception during simulation: {e}.")
        return np.inf  # Return a high value to indicate failure

if __name__ == "__main__":
    # Setup file logging
    logging.basicConfig(filename='differential_evolution_log.txt', level=logging.INFO, format='%(asctime)s - %(message)s')
    from random import random
    # Fixed parameters
    fixed_params = {
        'ifa_l': 14.0+random()*1.,  # Initial value
        'ifa_h': 7.5+random()*1.,
        'ifa_fp': 5.5+random()*1.,
        'ifa_w': 0.7+random()*0.1
    }

    # Define bounds for each variable you want to optimize
    variable_bounds = {
        'ifa_l': (10.0, 19.5),
        'ifa_h': (1.0, 14.0),
        'ifa_fp': (1.0, 10.0),
        'ifa_w': (0.4, 1.0)
    }

    # Variables to optimize
    variable_names = ['ifa_l', 'ifa_h', 'ifa_fp', 'ifa_w']
    bounds = [variable_bounds[var_name] for var_name in variable_names]

    # Run differential evolution
    result = differential_evolution(
        evaluation_fun,
        bounds,
        args=(variable_names, fixed_params),
        strategy='best1bin',
        maxiter=1000,
        popsize=15,
        tol=0.01,
        mutation=(0.5, 1),
        recombination=0.7,
        disp=True,
        polish=True
    )

    # Update fixed_params with the optimized values
    for i, var_name in enumerate(variable_names):
        fixed_params[var_name] = result.x[i]
        logging.info(f"Optimized {var_name}: {fixed_params[var_name]}")

    # Log the final result
    logging.info("Optimized parameters:")
    for var_name in variable_names:
        logging.info(f"{var_name}: {fixed_params[var_name]}")
    logging.info(f"Objective function value: {result.fun}")

    # Print the result to the console
    print("Optimization Result:")
    for var_name in variable_names:
        print(f"{var_name}: {fixed_params[var_name]}")
    print(f"Objective function value (S11 at center frequency): {result.fun}")
