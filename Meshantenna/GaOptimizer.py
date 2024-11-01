import random
import numpy as np
from GaAntenna import ga_simulation
from deap import base, creator, tools, algorithms
import logging
import os 
from time import time
import json
from scipy.ndimage import gaussian_filter  # Added for clustering

logpath = 'logs/ga_log800MHZ1800mhztest2.txt'
bestfitness = 0

# Define the shape of the 2D binary array
ARRAY_SHAPE = (20, 20)  # Ensure this shape fits your simulation requirements

# Number of elements in the 2D array
NUM_ELEMENTS = ARRAY_SHAPE[0] * ARRAY_SHAPE[1]

# Create the fitness and individual classes
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # Minimize S11
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# Remove the old biased attribute function
# import random

# def biased_attr_bool():
#     return random.choices([0, 1], weights=[0.6, 0.4])[0]

# Attribute generator for binary elements (0 or 1)
# toolbox.register("attr_bool", biased_attr_bool)

# New: Clustered antenna array generation function
def generate_clustered_array(array_size=ARRAY_SHAPE, mean=0.5, std_dev=0.3, correlation_length=2, threshold=0.501):
    # Generate random field
    random_field = np.random.normal(mean, std_dev, array_size)
    # Apply Gaussian filter to introduce spatial correlation
    smoothed_field = gaussian_filter(random_field, sigma=correlation_length)
    # Threshold the field to create the antenna array
    antenna_array = smoothed_field > threshold
    # Convert to int (0 or 1)
    antenna_array = antenna_array.astype(int)
    # Flatten the array to 1D
    flat_antenna_array = antenna_array.flatten()
    return flat_antenna_array.tolist()

# New: Custom individual initialization function
def init_clustered_individual():
    flat_antenna_array = generate_clustered_array()
    return creator.Individual(flat_antenna_array)

# Register the custom individual initialization function
toolbox.register("individual", init_clustered_individual)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# The rest of your code remains the same

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
        'plot': False,  # Set to True to plot results
        'showCad': False,
        'post_proc_only': False,
        'delete_simulation_files': True,
        'antenna_grid': array_2d,
        'randomhash': random.randint(0, 1000000),
        'frequencies': [0.83e9,1.8e9],
        'xmultiplier': 3,
        'ymultiplier': 3,
        'zmultiplier': 3,
        'lambdamultiplier': 2,
    }

    try:
        # Call the simulation function with the 2D array
        freq, s11_dB, Zin, P_in, hash_prefix = ga_simulation(params)

        # Extract S11 at center frequency
        fitness = 1.0
        s11_at_center, impedance, reactance = [], [], []

        for f in params["frequencies"]:
            idx = (np.abs(freq - f)).argmin()
            s11_value = s11_dB[idx]
            s11_at_center.append(s11_value)
            fitness *= abs(s11_value) if s11_value <= 0 else 1e-6  # Penalize positive S11 values

        fitness *= -1
        # bandwidth at xDB
        minS11 = -10

        impedance = np.real(Zin[idx])
        reactance = np.imag(Zin[idx])

        bandwidth = -1
        
        total_seconds = time() - starttime
        
        params_serializable = serialize_data(params)
        
        log_data = {
            "total_seconds": round(total_seconds, 2),
            "fitness": f"{fitness:.4f}",    
            "S11_at_cf": [f"{val:.3f}" for val in s11_at_center],
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
toolbox.register("mate", tools.cxUniform, indpb=0.1)  # Uniform crossover with low swap probability
toolbox.register("mutate", tools.mutFlipBit, indpb=0.02)  # Bit-flip mutation with lower probability
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
    population = toolbox.population(n=50)

    # Statistics to keep track of progress
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)

    # Hall of Fame to store the best individuals
    hof = tools.HallOfFame(1)

    # Genetic Algorithm parameters
    NGEN = 200  # Number of generations
    CXPB = 0.7  # Crossover probability
    MUTPB = 0.1  # Mutation probability

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
    logging.info(f"Best Fitness (S11): {hof[0].fitness.values[0]}")

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
