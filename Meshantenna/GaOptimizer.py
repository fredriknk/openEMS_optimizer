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
logpath = 'logs/ga_log800MHZ1800mhztest2.txt'

# Define the shape of the 2D binary array
ARRAY_SHAPE = (20, 20)

# Number of elements in the 2D array
NUM_ELEMENTS = ARRAY_SHAPE[0] * ARRAY_SHAPE[1]
NUM_PARAMS = 3  # std_dev, correlation_length, threshold
TOTAL_GENES = NUM_ELEMENTS + NUM_PARAMS

# Indices for the clustering parameters
std_dev_idx = NUM_ELEMENTS
correlation_length_idx = NUM_ELEMENTS + 1
threshold_idx = NUM_ELEMENTS + 2

# Create the fitness and individual classes
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # Minimize S11
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

def init_full_individual():
    # Initialize antenna elements
    antenna_elements = [random.choice([0, 1]) for _ in range(NUM_ELEMENTS)]
    
    # Define parameter ranges
    std_dev_range = (0.1, 1.0)
    correlation_length_range = (0.3, 5.0)
    threshold_range = (0.4, 0.59)
    
    # Initialize clustering parameters
    std_dev = random.uniform(*std_dev_range)
    correlation_length = random.uniform(*correlation_length_range)
    threshold = random.uniform(*threshold_range)
    
    # Combine antenna elements and clustering parameters
    individual_genes = antenna_elements + [std_dev, correlation_length, threshold]
    
    return creator.Individual(individual_genes)

def mutate_individual(individual):
    # Mutate antenna elements
    for i in range(NUM_ELEMENTS):
        if random.random() < 0.01:  # Mutation probability for antenna elements
            individual[i] = 1 - individual[i]  # Flip bit
    
    # Mutate clustering parameters with Gaussian mutation
    mutation_strength = 0.1
    std_dev_range = (0.1, 1.0)
    correlation_length_range = (1.0, 5.0)
    threshold_range = (0.3, 0.7)
    
    # Mutate std_dev
    if random.random() < 0.3:
        individual[std_dev_idx] += random.gauss(0, mutation_strength)
        individual[std_dev_idx] = min(max(individual[std_dev_idx], std_dev_range[0]), std_dev_range[1])
    # Mutate correlation_length
    if random.random() < 0.3:
        individual[correlation_length_idx] += random.gauss(0, mutation_strength)
        individual[correlation_length_idx] = min(max(individual[correlation_length_idx], correlation_length_range[0]), correlation_length_range[1])
    # Mutate threshold
    if random.random() < 0.3:
        individual[threshold_idx] += random.gauss(0, mutation_strength)
        individual[threshold_idx] = min(max(individual[threshold_idx], threshold_range[0]), threshold_range[1])
    
    return (individual,)

def mate_individuals(ind1, ind2):
    # Crossover antenna elements using uniform crossover
    for i in range(NUM_ELEMENTS):
        if random.random() < 0.5:
            ind1[i], ind2[i] = ind2[i], ind1[i]
    
    # Crossover clustering parameters using arithmetic crossover
    for idx in [std_dev_idx, correlation_length_idx, threshold_idx]:
        if random.random() < 0.5:
            ind1_val = ind1[idx]
            ind2_val = ind2[idx]
            ind1[idx] = (ind1_val + ind2_val) / 2.0
            ind2[idx] = (ind1_val + ind2_val) / 2.0
    return ind1, ind2

# Register the individual initialization and operators
toolbox.register("individual", init_full_individual)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", mate_individuals)
toolbox.register("mutate", mutate_individual)
toolbox.register("select", tools.selTournament, tournsize=3)

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
    
    # Extract antenna elements and clustering parameters
    flat_antenna_array = individual[:NUM_ELEMENTS]
    array_2d = np.array(flat_antenna_array).reshape(ARRAY_SHAPE)
    
    std_dev = individual[std_dev_idx]
    correlation_length = individual[correlation_length_idx]
    threshold = individual[threshold_idx]
    
    # Modify the antenna array using clustering parameters
    smoothed_array = gaussian_filter(array_2d.astype(float), sigma=correlation_length)
    array_2d = (smoothed_array > threshold).astype(int)
    flat_antenna_array = array_2d.flatten().tolist()
    
    # Update the individual's antenna elements
    individual[:NUM_ELEMENTS] = flat_antenna_array
    
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
            "S11_at_cf": [f"{val:.3f}" for val in s11_at_center],
            "std_dev": f"{std_dev:.3f}",
            "correlation_length": f"{correlation_length:.3f}",
            "threshold": f"{threshold:.3f}",
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
    CXPB = 0.7
    MUTPB = 0.2
    
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
    best_antenna_elements = best_individual[:NUM_ELEMENTS]
    best_array_2d = np.array(best_antenna_elements).reshape(ARRAY_SHAPE)
    std_dev = best_individual[std_dev_idx]
    correlation_length = best_individual[correlation_length_idx]
    threshold = best_individual[threshold_idx]
    
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
