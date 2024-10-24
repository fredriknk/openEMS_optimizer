from mifa_test import ifa_simulation
import logging
import numpy as np
from scipy.optimize import minimize, minimize_scalar

# Function to evaluate the antenna performance
def evaluation_fun(x, variable_names, fixed_params):
    # Update the variables to be optimized
    params = fixed_params.copy()
    
    if len(variable_names) == 1:
        # Ensure x is a scalar
        if isinstance(x, (list, np.ndarray)):
            params[variable_names[0]] = x[0]
        else:
            params[variable_names[0]] = x
    else:
        # Treat x as a vector
        for i, var_name in enumerate(variable_names):
            params[var_name] = x[i]

    # Debugging: Verify parameter types
    for key, value in params.items():
        print(f"{key}: {value}, type: {type(value)}")

    # Extract parameters
    ifa_h = params['ifa_h']
    ifa_l = params['ifa_l']
    ifa_fp = params['ifa_fp']
    ifa_w1 = params['ifa_w1']
    ifa_w2 = params['ifa_w2']
    ifa_wf = params['ifa_wf']

    # Constants for the simulation
    Sim_CSX = 'IFA.xml'
    showCad = False
    post_proc_only = False
    unit = 1e-3
    substrate_width = 21
    substrate_length = 40
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
        freq, s11_dB, Zin, P_in,hash_prefix = ifa_simulation(
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
            ifa_w1=ifa_w1,
            ifa_w2=ifa_w2,
            ifa_wf=ifa_wf,
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
        #round value to 2 decimal places
        
        f_res_ind = np.argmin(s11_dB)
        f_res = freq[f_res_ind]
        s_11_min = s11_dB[f_res_ind]
        #bandwidth at xDB
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
        #specify the bandwidth:
        if first_crossing is not -1 or last_crossing is not -1:
            bandwidth = last_crossing - first_crossing  
        
        # Log parameters and objective function values
        res = 4
        logging.info(f"ifa_l: {ifa_l:.3f}, ifa_h: {ifa_h:.3f}, ifa_fp: {ifa_fp:.3f}, ifa_w1: {ifa_w1:.3f},ifa_w1: {ifa_w2:.3f},ifa_wf: {ifa_wf:.3f}, S11 at cf: {s11_value:.4f}, Imp:{impedance:.1f}R {reactance:.1f}z, Res f:{f_res/1e9:.3f} GHz, S11 at res:{s_11_min:.3f}, BW1:{first_crossing/1e9:.2f} GHz, BW2:{last_crossing/1e9:.2f} GHz, BW = {bandwidth/1e6:.1f} MHz - id {hash_prefix}")
        # Return the magnitude (since we want to minimize it)
        return s11_value

    except Exception as e:
        logging.info(f"Exception during simulation: {e}.")
        return np.inf  # Return a high value to indicate failure


if __name__ == "__main__":
    # Setup file logging
    logging.basicConfig(filename='logs\\gradient_bisect_log.txt', level=logging.INFO, format='%(asctime)s - %(message)s')

    # Fixed parameters
    fixed_params = {
        'ifa_l': 23,  # Initial value
        'ifa_h': 5.5,  
        'ifa_fp': 5,
        'ifa_w1': 1.5,
        'ifa_w2': 0.5,
        'ifa_wf': 0.5,
    }

    # Define bounds for each variable you want to optimize
    variable_bounds = {
        'ifa_l': (16, 31),
        'ifa_h': (5., 16.0),
        'ifa_fp': (1.0, 6),
        'ifa_w1': (0.35, 3),
        'ifa_w2': (0.4, 1.0),
        'ifa_wf': (0.4, 1.5)
    }
    
    # Choose variables to optimize
    variable_names = ['ifa_l', 'ifa_w1', 'ifa_l','ifa_w2','ifa_wf','ifa_fp','ifa_h']  # List variables you want to optimize
    bcplist = variable_names

    for var_name in bcplist:
        logging.info(f"Optimizing {var_name} within bounds: {variable_bounds[var_name]}")
        variable_names = [var_name]
        bounds = [variable_bounds[var_name]]
        fixed_params_copy = fixed_params.copy()
        print(f"bounds: {bounds}")
        
        result = minimize_scalar(
            evaluation_fun,
            bounds=(bounds[0][0], bounds[0][1]),
            method='bounded',
            args=(variable_names, fixed_params_copy),
            options={'xatol': (bounds[0][1] - bounds[0][0]) / 100}  # Stop when the interval is smaller than 1/100th the bounds
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