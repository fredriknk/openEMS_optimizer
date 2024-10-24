from mifa_test import ifa_simulation
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
    min_size = 0.3  # Minimum automesh size
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
        s11_value = round(s11_value, 2)
        
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
        logging.info(f"ifa_l: {ifa_l:.3f}, ifa_h: {ifa_h:.3f}, ifa_fp: {ifa_fp:.3f}, ifa_w1: {ifa_w1:.3f},ifa_w1: {ifa_w2:.3f},ifa_wf: {ifa_wf:.3f}, S11 at cf: {s11_value}, Imp:{impedance:.1f}R {reactance:.1f}z, Res f:{f_res/1e9:.3f} GHz, S11 at res:{s_11_min:.3f}, BW1:{first_crossing/1e9:.2f} GHz, BW2:{last_crossing/1e9:.2f} GHz, BW = {bandwidth/1e6:.1f} MHz - id {hash_prefix}")
        # Return the magnitude (since we want to minimize it)
        return s11_value

    except Exception as e:
        logging.info(f"Exception during simulation: {e}.")
        return np.inf  # Return a high value to indicate failure


if __name__ == "__main__":
    # Setup file logging
    logging.basicConfig(filename='logs\\optimization_log.txt', level=logging.INFO, format='%(asctime)s - %(message)s')

    # Fixed parameters
    fixed_params = {
        'ifa_l': 20,  # Initial value
        'ifa_h': 7.5,  
        'ifa_fp': 7,
        'ifa_w1': 0.75,
        'ifa_w2': 0.75,
        'ifa_wf': 0.75,
    }

    # Define bounds for each variable you want to optimize
    variable_bounds = {
        'ifa_l': (16, 31),
        'ifa_h': (5., 16.0),
        'ifa_fp': (1.0, 7),
        'ifa_w1': (0.4, 1.0),
        'ifa_w2': (0.4, 1.0),
        'ifa_wf': (0.6, 1.5),
    }

    # Choose variables to optimize
    variable_names = ['ifa_l', 'ifa_h', 'ifa_fp', 'ifa_wf','ifa_w2','ifa_w1']  # List variables you want to optimize
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

        # Use the best initial value as x0 for the optimizer
        x0 = [best_initial_value]
        logging.info(f"Optimizing {var_name} within bounds: {variable_bounds[var_name]} starting from x0: {x0}")
        result = minimize(
            evaluation_fun,
            x0,
            args=(variable_names, fixed_params),
            method='L-BFGS-B',
            bounds=bounds,
            options={'disp': True, 'maxiter': 10,'maxfun': 5, "eps":0.05,'gtol': 1e-2,'ftol': 1e-2}
        )

        fixed_params[var_name] = result.x[0]
        logging.info(f"Optimized {var_name}: {fixed_params[var_name]}, Objective function value: {result.fun}")
        
    # Choose variables to optimize
    variable_names = ['ifa_fp', 'ifa_l', 'ifa_h', 'ifa_fp', 'ifa_wf','ifa_w2','ifa_w1']  # List variables you want to optimize further
    
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
    
        # Use the best initial value as x0 for the optimizer
        x0 = fixed_params[var_name]
        logging.info(f"Optimizing {var_name} within bounds: {variable_bounds[var_name]} starting from x0: {x0}")
        result = minimize(
            evaluation_fun,
            x0,
            args=(variable_names, fixed_params),
            method='L-BFGS-B',
            bounds=bounds,
            options={'disp': True, 'maxiter': 10,'maxfun': 5, "eps":1e-4}#,'gtol': 1e-5,'ftol': 1e-4}
        )

        fixed_params[var_name] = result.x[0]
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
