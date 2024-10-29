from mifa_groundplane import ifa_simulation
import logging
import numpy as np
import os
import pickle
from scipy.optimize import differential_evolution


# Function to save the state of the differential evolution
def save_state(filename, state):
    with open(filename, 'wb') as f:
        pickle.dump(state, f)
    logging.info(f"State saved to {filename}")


# Function to load the state of the differential evolution
def load_state(filename):
    if os.path.exists(filename):
        with open(filename, 'rb') as f:
            state = pickle.load(f)
        logging.info(f"State loaded from {filename}")
        return state
    return None


# Callback function to save state periodically
def save_callback(xk, convergence):
    global result_state, save_filename
    # Save the current state of the optimization
    result_state = {'xk': xk, 'convergence': convergence, 'population': optimizer.population}
    save_state(save_filename, result_state)
    return False  # Returning False allows the optimization to continue


def evaluation_fun(x, variable_names, fixed_params):
    import logging
    import os
    import numpy as np
    import traceback
    from time import time
    starttime = time()
    # Ensure the logs directory exists
    os.makedirs('logs', exist_ok=True)

    # Get the root logger
    logger = logging.getLogger()
    if not logger.hasHandlers():
        # Configure logging only if not already configured
        logging.basicConfig(
            filename='logs/diffevolution_log.txt',
            level=logging.INFO,
            format='%(asctime)s - %(message)s',
            filemode='a'  # Append mode
        )

    # Initialize params by copying fixed_params
    params = fixed_params.copy()

    # Update the variables to be optimized
    if len(variable_names) == 1:
        # Ensure x is a scalar
        params[variable_names[0]] = x if np.isscalar(x) else x[0]
    else:
        # Treat x as a vector
        for i, var_name in enumerate(variable_names):
            params[var_name] = x[i]

    # Optional: Log parameter types for debugging
    for key, value in params.items():
        logger.debug(f"{key}: {value}, type: {type(value)}")

    # Extract parameters
    ifa_h = params['ifa_h']
    ifa_l = params['ifa_l']
    ifa_fp = params['ifa_fp']
    ifa_w1 = params['ifa_w1']
    ifa_w2 = params['ifa_w2']
    ifa_wf = params['ifa_wf']
    center_freq = 2.45e9  # Center frequency for S11 evaluation
    # Constants for the simulation (as per your original code)
    # ...

    try:
        freq, s11_dB, Zin, P_in, hash_prefix = ifa_simulation(
            Sim_CSX='IFA.xml',
            showCad=False,
            post_proc_only=False,
            unit=1e-3,
            substrate_width=21,
            substrate_length=83,
            substrate_thickness=1.5,
            gndplane_position=0,
            substrate_cells=4,
            ifa_h=ifa_h,
            ifa_l=ifa_l,
            ifa_fp=ifa_fp,
            ifa_w1=ifa_w1,
            ifa_w2=ifa_w2,
            ifa_wf=ifa_wf,
            ifa_e=0.5,
            substrate_epsR=4.5,
            feed_R=50,
            min_freq=2.4e9,
            center_freq=center_freq,
            max_freq=2.5e9,
            plot=False,
            min_size=0.3
        )

        # Get the S11 value at the center frequency
        idx = np.argmin(np.abs(freq - center_freq))
        s11_value = s11_dB[idx]
        # round value to 2 decimal places

        f_res_ind = np.argmin(s11_dB)
        f_res = freq[f_res_ind]
        s_11_min = s11_dB[f_res_ind]
        # bandwidth at xDB
        minS11 = -6

        first_crossing = -1
        last_crossing = -1
        if f_res_ind != 0 and s_11_min < minS11:
            for i in range(0, f_res_ind):
                if s11_dB[i] < minS11:
                    first_crossing = freq[i]
                    break

            # Find the last frequency crossing over -10 dB after f_res_ind

            for i in range(f_res_ind, len(s11_dB)):
                if s11_dB[i] > minS11:
                    last_crossing = freq[i]
                    break

        impedance = np.real(Zin[idx])
        reactance = np.imag(Zin[idx])

        bandwidth = -1
        # specify the bandwidth:
        if first_crossing is not -1 or last_crossing is not -1:
            bandwidth = last_crossing - first_crossing
        total_seconds = time() - starttime
        # Log parameters and objective function values
        log_message = (
            f"total seconds: {total_seconds:.2f}, "
            f"ifa_l: {ifa_l:.3f}, ifa_h: {ifa_h:.3f}, ifa_fp: {ifa_fp:.3f}, "
            f"ifa_w1: {ifa_w1:.3f}, ifa_w2: {ifa_w2:.3f}, ifa_wf: {ifa_wf:.3f}, "
            f"S11 at cf: {s11_value:.4f}, Imp: {impedance:.1f}R {reactance:.1f}z, "
            f"Res f: {f_res / 1e9:.3f} GHz, S11 at res: {s_11_min:.3f}, "
            f"BW1: {first_crossing / 1e9:.2f} GHz, BW2: {last_crossing / 1e9:.2f} GHz, "
            f"BW = {bandwidth / 1e6:.1f} MHz - id {hash_prefix}"
        )
        logger.info(log_message)

        # Return the objective function value (since we want to minimize it)
        return s11_value

    except Exception as e:
        logger.error(f"Exception in evaluation_fun: {e}")
        logger.error(traceback.format_exc())
        return np.inf


if __name__ == "__main__":
    # Ensure the logs directory exists
    os.makedirs('logs', exist_ok=True)
    os.makedirs("savefiles", exist_ok=True)
    save_filename = 'savefiles/diffevolution_state.pkl'

    # Configure logging once in the main process
    logging.basicConfig(
        filename='logs/diffevolution_log.txt',
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        filemode='a'  # Append mode
    )

    # Fixed parameters
    fixed_params = {
        'ifa_l': 23,  # Initial value
        'ifa_h': 5.5,
        'ifa_fp': 5,
        'ifa_w1': 1,
        'ifa_w2': 1.,
        'ifa_wf': 1,
    }

    # Define bounds for each variable you want to optimize
    variable_bounds = {
        'ifa_l': (19, 24),
        'ifa_h': (5., 7.0),
        'ifa_fp': (4.0, 6),
        'ifa_w1': (0.5, 3),
        'ifa_w2': (0.4, 1.5),
        'ifa_wf': (0.4, 1.5)
    }

    # Variables to optimize
    variable_names = [ 'ifa_w1', 'ifa_l','ifa_w2','ifa_wf','ifa_fp','ifa_h']  # List variables you want to optimize
    bounds = [variable_bounds[var_name] for var_name in variable_names]

    logging.info(f"start diff evolution, bounds: {bounds}, fixed_params: {fixed_params}")

    # Check for saved state
    result_state = load_state(save_filename)
    init_pop = None

    # If a saved state exists, extract the initial population
    if result_state:
        init_pop = result_state.get('population')

    # Run differential evolution with the initial population if available
    optimizer = differential_evolution(
        evaluation_fun,
        bounds,
        args=(variable_names, fixed_params),
        strategy='best1bin',
        maxiter=1000,
        popsize=5,
        tol=0.01,
        mutation=(0.5, 1),
        recombination=0.7,
        disp=True,
        polish=True,
        workers=4,
        init=init_pop if init_pop is not None else 'random'
    )

    # After completion, save the final state
    result_state = {'xk': optimizer.x, 'convergence': optimizer.fun, 'population': optimizer.population}
    save_state(save_filename, result_state)

    # Update fixed_params with the optimized values
    for i, var_name in enumerate(variable_names):
        fixed_params[var_name] = optimizer.x[i]
        logging.info(f"Optimized {var_name}: {fixed_params[var_name]}")

    # Log the final result
    logging.info("Optimized parameters:")
    for var_name in variable_names:
        logging.info(f"{var_name}: {fixed_params[var_name]}")
    logging.info(f"Objective function value: {optimizer.fun}")

    # Print the result to the console
    print("Optimization Result:")
    for var_name in variable_names:
        print(f"{var_name}: {fixed_params[var_name]}")
    print(f"Objective function value (S11 at center frequency): {optimizer.fun}")
