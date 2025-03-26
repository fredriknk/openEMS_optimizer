from mifa import ifa_simulation
import logging
import numpy as np
import os
import pickle
from scipy.optimize import differential_evolution
from datetime import datetime

time_now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S") 

logname = f"logs/diffevolution_log_{time_now}.txt"
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

def reflection_percent(s11_db):
    """
    Calculate the percentage of incident power that is reflected.

    :param s11_db: (float) S11 in dB
    :return: (float) Reflected power in percent
    """
    gamma_magnitude = 10 ** (s11_db / 20.0)   # Reflection coefficient |Γ|
    return 1- (gamma_magnitude ** 2)     # % reflected = 100 * |Γ|^2

def evaluation_fun(x, variable_names, fixed_params,logname):
    import logging
    import os
    import numpy as np
    import traceback
    from time import time
    starttime = time()
    # Ensure the logs directory exists
    os.makedirs('logs', exist_ok=True)
    
    Sim_CSX = 'IFA.xml'
    showCad = False
    post_proc_only = False
    unit = 1e-3
    substrate_width = 25.5
    substrate_length = 108
    substrate_thickness = 1.5
    gndplane_position = -0.36
    substrate_cells = 4
    #ifa_h = 16
    #ifa_l = 200
    #ifa_w1 = 2
    #ifa_w2 = 0.5
    #ifa_wf = 0.5
    #ifa_fp = 4.5
    ifa_e = 0.5
    mifa_meander=2.2
    mifa_tipdistance=2
    mifa_meander_edge_distance=mifa_tipdistance
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 0.791e9
    center_freq = 0.826e9
    max_freq = 0.862e9
    fc = 0.83e9-0.7e9
    min_size = 0.3 # minimum automesh sizee
    max_size=2 #maximum automesh size
    override_min_global_grid = None #none if not override
    max_timesteps = 1000000
    plot = False
    cleanup = False
    
    # Get the root logger
    logger = logging.getLogger()
    if not logger.hasHandlers():
        # Configure logging only if not already configured
        logging.basicConfig(
            filename=logname,
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
    ifa_l = params['ifa_l']
    ifa_h = params['ifa_h']
    ifa_fp = params['ifa_fp']
    ifa_w1 = params['ifa_w1']
    ifa_w2 = params['ifa_w2']
    ifa_wf = params['ifa_wf']
    
    # Constants for the simulation (as per your original code)
    # ...

    try:
        freq, s11_dB, Zin, P_in, hash_prefix,efficiency = ifa_simulation(Sim_CSX=Sim_CSX,
                                                showCad=showCad,
                                                post_proc_only=post_proc_only,
                                                unit=post_proc_only,
                                                substrate_width=substrate_width,
                                                substrate_length=substrate_length,
                                                substrate_thickness=substrate_thickness,
                                                gndplane_position=gndplane_position,
                                                substrate_cells=substrate_cells,
                                                ifa_h=ifa_h,
                                                ifa_l=ifa_l,
                                                ifa_w1=ifa_w1,
                                                ifa_w2=ifa_w2,
                                                ifa_wf=ifa_wf,
                                                ifa_fp=ifa_fp,
                                                ifa_e=ifa_e,
                                                mifa_meander=mifa_meander,
                                                mifa_tipdistance=mifa_tipdistance,
                                                mifa_meander_edge_distance=mifa_meander_edge_distance,
                                                substrate_epsR=substrate_epsR,
                                                feed_R=feed_R,
                                                min_freq=min_freq,
                                                center_freq=center_freq,
                                                fc=fc,
                                                max_freq=max_freq,
                                                min_size=min_size,
                                                max_size=max_size,
                                                override_min_global_grid=override_min_global_grid,
                                                max_timesteps=max_timesteps,
                                                plot=plot,
                                                delete_simulation_files=cleanup)
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
            f"Efficiency: {efficiency:.3f}, "
            f"Comps11 at cf: {reflection_percent(s11_value)*(efficiency/100):.4f}, "
            f"S11 at cf: {s11_value:.4f}, Imp: {impedance:.1f}R {reactance:.1f}z, "
            f"ifa_l: {ifa_l:.3f}, ifa_h: {ifa_h:.3f}, ifa_fp: {ifa_fp:.3f}, "
            f"ifa_w1: {ifa_w1:.3f}, ifa_w2: {ifa_w2:.3f}, ifa_wf: {ifa_wf:.3f}, "
            f"mifa_meander_edge_distance: {mifa_meander_edge_distance:.3f}, "
            f"Res f: {f_res / 1e9:.3f} GHz, S11 at res: {s_11_min:.3f}, "
            f"BW1: {first_crossing / 1e9:.2f} GHz, BW2: {last_crossing / 1e9:.2f} GHz, "
            f"BW = {bandwidth / 1e6:.1f} MHz - id {hash_prefix}"
        )
        logger.info(log_message)

        # Return the objective function value (since we want to minimize it)
        return reflection_percent(s11_value)*(efficiency/100)

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
        filename=logname,
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        filemode='a'  # Append mode
    )

    # Fixed parameters
    fixed_params = {
        "ifa_l" : 85,
        "ifa_h" : 13,
        "ifa_w1" : 0.6,
        "ifa_w2" : 0.6,
        "ifa_wf" : 0.6,
        "ifa_fp" : 2.914,
    }

    # Define bounds for each variable you want to optimize
    variable_bounds = {
        'ifa_l': (60, 94),
        'ifa_h': (12, 16),
        'ifa_fp': (2.5, 6),
        'ifa_w1': (0.4, 1.5),
        'ifa_w2': (0.4, 1.0),
        'ifa_wf': (0.4, 1.0),
    }

    # Variables to optimize from variable_bounds
    variable_names = list(variable_bounds.keys())   
    #variable_names = [ 'ifa_l',"ifa_h"]#,'ifa_fp', 'ifa_w1','ifa_w2','ifa_wf']  # List variables you want to optimize
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
        args=(variable_names, fixed_params,logname),
        strategy='best1bin',
        maxiter=1000,
        popsize=7,
        tol=0.1,
        mutation=(0.5, 1),
        recombination=0.7,
        disp=True,
        polish=True,
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
