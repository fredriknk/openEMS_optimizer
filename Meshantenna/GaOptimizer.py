import random
import numpy as np
from GaAntenna import ga_simulation
from deap import base, creator, tools, algorithms
import logging
import os 
from time import time
import json


logpath = 'logs/ga_log.txt'
bestfitness = 0

# Define the shape of the 2D binary array
ARRAY_SHAPE = (40, 40)  # Example shape, adjust as needed

# Number of elements in the 2D array
NUM_ELEMENTS = ARRAY_SHAPE[0] * ARRAY_SHAPE[1]

# Create the fitness and individual classes
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # Minimize S11
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

import random

def biased_attr_bool():
    return random.choices([0, 1], weights=[0.6, 0.4])[0]

# Attribute generator for binary elements (0 or 1)
toolbox.register("attr_bool", biased_attr_bool)

# Structure initializer for individuals (flattened binary arrays)
toolbox.register(
    "individual",
    tools.initRepeat,
    creator.Individual,
    toolbox.attr_bool,
    n=NUM_ELEMENTS,
)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)



def serialize_data(data):
    if isinstance(data, np.ndarray):
        return data.tolist()  # Convert ndarray to list
    elif isinstance(data, dict):
        return {key: serialize_data(value) for key, value in data.items()}  # Recursively handle dict
    elif isinstance(data, list):
        return [serialize_data(item) for item in data]  # Recursively handle list
    else:
        return data  # Base case: return the item itself if it's already JSON serializable


def evaluate(individual):
    starttime = time()
     # Get the root logger
    logger = logging.getLogger()
    if not logger.hasHandlers():
        # Configure logging only if not already configured
        logging.basicConfig(
            filename=logpath,
            level=logging.INFO,
            format='%(asctime)s - %(message)s',
            filemode='a'  # Append mode
    )
    # Reshape the individual back to the original 2D array shape
    array_2d = np.array(individual).reshape(ARRAY_SHAPE)

    # Fixed parameters
    params = {
        'Sim_CSX' : 'IFA.xml',
        'unit': 1e-3,
        'substrate_width': 21,
        'substrate_length': 20,
        'substrate_thickness': 1.5,
        'substrate_epsR': 4.5,
        'gndplane_position': 0,
        'substrate_cells': 4,
        'ant_h': 14,
        'ant_l': 20,
        'ant_fp': 5,
        'ant_e': 0.5,
        'feed_R': 50,
        'min_freq': 2.4e9,
        'center_freq': 2.45e9,
        'max_freq': 2.5e9,
        'fc': 1.0e9,
        'max_timesteps': 12000,
        'override_min_global_grid': None,
        'plot': False,  # Set to True to plot results
        'showCad': False,
        'post_proc_only': False,
        'delete_simulation_files': True,
        'antenna_grid': array_2d,
        'randomhash': random.randint(0, 1000000)
    }

    try:
        # Call the simulation function with the 2D array
        freq, s11_dB, Zin, P_in, hash_prefix = ga_simulation(params)

        # Extract S11 at center frequency
        center_freq = params['center_freq']
        idx = (np.abs(freq - center_freq)).argmin()
        s11_at_center = s11_dB[idx]
        fitness = s11_at_center
        
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
        
        params_serializable = serialize_data(params)
        
        log_data = {
            "total_seconds": round(total_seconds, 2),
            "S11_at_cf": round(s11_value, 4),
            "impedance": f"{impedance:.1f}R",
            "reactance": f"{reactance:.1f}z",
            "resonant_frequency": f"{f_res / 1e9:.3f} GHz",
            "S11_at_resonant_frequency": round(s_11_min, 3),
            "bandwidth_lower": round(first_crossing / 1e9, 2),
            "bandwidth_upper": round(last_crossing / 1e9, 2),
            "bandwidth": f"{bandwidth / 1e6:.1f} MHz",
            "params": params_serializable  # Convert array to list
        }
        log_message = json.dumps(log_data, separators=(',', ':'))
        # Log the results as a single line of JSON
        logging.info(log_message)

    except Exception as e:
        
        # In case of simulation failure, assign a high fitness value
        print(f"Simulation failed for individual {individual}: {e}")
        logging.error(f"Simulation failed for individual {individual}: {e}")
        fitness = 100.0  # Assign a large penalty

    return (fitness,)

toolbox.register("evaluate", evaluate)

# Genetic Algorithm Operators
toolbox.register("mate", tools.cxTwoPoint)  # Two-point crossover
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)  # Bit-flip mutation
toolbox.register("select", tools.selTournament, tournsize=3)

def optimize_antenna():
    os.makedirs('logs', exist_ok=True)
    # Configure logging once in the main process
    logging.basicConfig(
        filename=logpath,
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        filemode='a'  # Append mode
    )
    # Initialize population
    population = toolbox.population(n=30)

    # Statistics to keep track of progress
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)

    # Hall of Fame to store the best individuals
    hof = tools.HallOfFame(1)

    # Genetic Algorithm parameters
    NGEN = 100  # Number of generations
    CXPB = 0.5  # Crossover probability
    MUTPB = 0.2  # Mutation probability

    # Run the Genetic Algorithm
    algorithms.eaSimple(
        population,
        toolbox,
        cxpb=CXPB,
        mutpb=MUTPB,
        ngen=NGEN,
        stats=stats,
        halloffame=hof,
        verbose=True,
    )

    # Best individual
    best_individual = hof[0]
    best_array_2d = np.array(best_individual).reshape(ARRAY_SHAPE)
    logging.info("Best 2D Binary Array:")
    logging.info(best_array_2d)
    logging.info("Best Fitness (S11):", hof[0].fitness.values[0])

    return best_array_2d

if __name__ == "__main__":
    best_array_2d = optimize_antenna()
    # Run a final simulation with the best 2D array
    params = {
        'Sim_CSX' : 'IFA.xml',
        'unit': 1e-3,
        'substrate_width': 21,
        'substrate_length': 20,
        'substrate_thickness': 1.5,
        'substrate_epsR': 4.5,
        'gndplane_position': 0,
        'substrate_cells': 4,
        'ant_h': 14,
        'ant_l': 20,
        'ant_fp': 5,
        'ant_e': 0.5,
        'feed_R': 50,
        'min_freq': 2.4e9,
        'center_freq': 2.45e9,
        'max_freq': 2.5e9,
        'fc': 1.0e9,
        'max_timesteps': 600000,
        'override_min_global_grid': None,
        'plot': False,  # Set to True to plot results
        'showCad': True,
        'post_proc_only': False,
        'delete_simulation_files': True,
        'antenna_grid': best_array_2d
    }
    freq, s11_dB, Zin, P_in, hash_prefix = ga_simulation(params)
    # Now you can plot or analyze the results as needed
