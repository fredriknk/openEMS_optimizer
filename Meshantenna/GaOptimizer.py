import random
import numpy as np
from GaAntenna import ga_simulation
from deap import base, creator, tools, algorithms
import logging
import os 
from time import time
import json
from scipy.ndimage import gaussian_filter
import math as m
logpath = 'logs/ga_log800MHZ1800mhztest3.txt'

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
    r = random.random()
    ri = 1-r
    return random.choices([0, 1], weights=[r, ri])[0]

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
        return data.tolist()
    elif isinstance(data, dict):
        return {key: serialize_data(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [serialize_data(item) for item in data]
    else:
        return data

def evaluate(individual):
    starttime = time()
    logger = logging.getLogger()
    if not logger.hasHandlers():
        logging.basicConfig(
            filename=logpath,
            level=logging.INFO,
            format='%(asctime)s - %(message)s',
            filemode='a'
        )
    
    # Reshape the individual back to the original 2D array shape
    array_2d = np.array(individual).reshape(ARRAY_SHAPE)
    
    # Fixed parameters
    params = {
        'Sim_CSX' : 'IFA.xml',
        'unit': 1e-3,
        'substrate_width': 25.5,
        'substrate_length': 107,
        'substrate_thickness': 1.5,
        'substrate_epsR': 4.5,
        'gndplane_position': 0,
        'substrate_cells': 4,
        'ant_h': 14,
        'ant_l': 24.5,
        'ant_fp': 5,
        'ant_e': 0.5,
        'feed_R': 50,
        'min_freq': 0.83e9,
        'center_freq': 1.5e9,
        'max_freq': 1.8e9,
        'fc': 0.8e9,
        'max_timesteps': 20000,
        'override_min_global_grid': None,
        'plot': False,
        'showCad': False,
        'post_proc_only': False,
        'delete_simulation_files': True,
        'antenna_grid': array_2d,
        'randomhash': random.randint(0, 1000000),
        'frequencies': [0.83e9, 1.8e9],
        'numThreads': 4,
        'xmultiplier': 3,
        'ymultiplier': 3,
        'zmultiplier': 3,
        'lambdamultiplier': 2,
    }
    
    try:
        # Call the simulation function with the modified 2D array
        freq, s11_dB, Zin, P_in, hash_prefix = ga_simulation(params)
        
        # Extract S11 at center frequency
        fitness = 1.0
        s11_at_center = []
        
        for f in params["frequencies"]:
            idx = (np.abs(freq - f)).argmin()
            s11_value = s11_dB[idx]
            s11_at_center.append(s11_value)
            fitness *= abs(s11_value) if s11_value <= 0 else 1e-6
        fitness = m.sqrt(fitness)   
        fitness *= -1
        
        total_seconds = time() - starttime
        
        params_serializable = serialize_data(params)
        
        log_data = {
            "total_seconds": round(total_seconds, 2),
            "fitness": f"{fitness:.5f}",
            "S11_at_f": [f"{val:.3f}" for val in s11_at_center],
            "params": params_serializable
        }
        log_message = json.dumps(log_data, separators=(',', ':'))
        logging.info(log_message)
    
    except Exception as e:
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
    logging.basicConfig(
        filename=logpath,
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        filemode='a'
    )
    # Initialize population
    population = toolbox.population(n=50)
    
    # Statistics to keep track of progress
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)
    
    # Hall of Fame to store the best individuals
    hof = tools.HallOfFame(1)
    
    # Genetic Algorithm parameters
    NGEN = 200
    CXPB = 0.8
    MUTPB = 0.05
    
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
    logging.info(f"Best Clustering Parameters: std_dev={std_dev}, correlation_length={correlation_length}, threshold={threshold}")
    logging.info(f"Best Fitness (S11): {hof[0].fitness.values[0]}")
    
    return best_array_2d

if __name__ == "__main__":
    array_2d = optimize_antenna()
    # Run a final simulation with the best 2D array
    params = {
        'Sim_CSX' : 'IFA.xml',
        'unit': 1e-3,
        'substrate_width': 25,
        'substrate_length': 80,
        'substrate_thickness': 1.5,
        'substrate_epsR': 4.5,
        'gndplane_position': 0,
        'substrate_cells': 4,
        'ant_h': 14,
        'ant_l': 20,
        'ant_fp': 5,
        'ant_e': 0.5,
        'feed_R': 50,
        'min_freq': 0.83e9,
        'center_freq': 1.5e9,
        'max_freq': 1.8e9,
        'fc': 0.8e9,
        'max_timesteps': 20000,
        'override_min_global_grid': None,
        'plot': False,
        'showCad': False,
        'post_proc_only': False,
        'delete_simulation_files': True,
        'antenna_grid': array_2d,
        'randomhash': random.randint(0, 1000000),
        'frequencies': [0.83e9, 1.8e9],
        'xmultiplier': 3,
        'ymultiplier': 3,
        'zmultiplier': 3,
        'lambdamultiplier': 2,
    }
    freq, s11_dB, Zin, P_in, hash_prefix = ga_simulation(params)
    # Now you can plot or analyze the results as needed
