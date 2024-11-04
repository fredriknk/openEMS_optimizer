import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from deap import base, creator, tools, algorithms
from scipy.ndimage import convolve

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
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))  # Minimize some objective
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

std_dev_range = (0.2, 0.25)
correlation_length_range = (0.4, 3.0)
threshold_range = (0.505, 0.52)

# Define the shape of the 2D binary array
ARRAY_SHAPE = (40, 40)
NUM_ELEMENTS = ARRAY_SHAPE[0] * ARRAY_SHAPE[1]

# Genetic Algorithm parameters 
NGEN = 5        # Number of generations
POP_SIZE = 10   # Population size
CXPB = 0.7      # Crossover probability
MUTPB = 0.2     # Mutation probability
BLOCK_SIZE = (5, 5)  # Size of the block for crossover and mutation


def generate_clustered_array(array_size=(20, 20), p_zero=0.5, seed=None):
    """
    Generates a 2D grid of 0s and 1s where each 1 has at least one orthogonal neighboring 1.
    
    Args:
        array_size (tuple of int, optional): Size of the grid as (rows, cols). Defaults to (40, 40).
        p_zero (float, optional): Probability of a cell being 0. Must be between 0 and 1. Defaults to 0.5.
        seed (int, optional): Seed for the random number generator for reproducibility. Defaults to None.
    
    Returns:
        list of int: Flattened 1D list representing the 2D grid with no isolated 1s.
    
    Raises:
        ValueError: If p_zero is not between 0 and 1.
    """
    rows, cols = array_size
    p_one = 1 - p_zero  # Probability of a cell being 1
    
    # Validate probabilities
    if not (0 <= p_zero <= 1) or not (0 <= p_one <= 1):
        raise ValueError("p_zero and p_one must be between 0 and 1.")
    
    # Seed the random number generators for reproducibility
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Step 1: Generate the initial grid with 0s and 1s based on probabilities
    grid = np.array([
        [random.choices([0, 1], weights=[p_zero, p_one])[0] for _ in range(cols)]
        for _ in range(rows)
    ])
    
    # Step 2: Define the convolution kernel to count orthogonal neighbors only
    kernel = np.array([
        [1, 1, 1],
        [1, 0, 1],
        [1, 1, 1]
    ])
    
    # Step 3: Convolve the grid with the kernel to count the number of orthogonal 1s around each cell
    neighbor_count = convolve(grid, kernel, mode='constant', cval=0)
    
    # Step 4: Identify isolated 1s (1s with no orthogonal neighboring 1s)
    isolated = (grid == 1) & (neighbor_count == 0)
    
    # Step 5: Set isolated 1s to 0
    grid[isolated] = 0
    
    # Optional: Repeat the process to remove any new isolated 1s resulting from previous removals
    # This ensures that all 1s have at least one orthogonal neighbor
    # Uncomment the following lines if multiple passes are necessary
    while True:
        neighbor_count = convolve(grid, kernel, mode='constant', cval=0)
        new_isolated = (grid == 1) & (neighbor_count == 0)
        if not new_isolated.any():
            break
        grid[new_isolated] = 0
    
    # Flatten the 2D grid to a 1D list and return
    flat_antenna_array = grid.flatten().tolist()
    return flat_antenna_array

def init_full_individual():
    """
    Initialize an individual with a clustered antenna array.
    """
    antenna_elements = generate_clustered_array(
        array_size=ARRAY_SHAPE
    )
    return creator.Individual(antenna_elements)

toolbox.register("individual", init_full_individual)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

from scipy.ndimage import gaussian_filter
import numpy as np

toolbox.register("mate", tools.cxTwoPoint)  # Two-point crossover
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)  # Bit-flip mutation
toolbox.register("select", tools.selTournament, tournsize=3)

# Evaluation function (dummy example)
def eval_individual(individual):
    # Dummy fitness function: count number of ones in antenna elements
    antenna_elements = individual[:NUM_ELEMENTS]
    fitness = antenna_elements.count(1)
    return (fitness,)

toolbox.register("evaluate", eval_individual)

# Genetic Algorithm parameters 
NGEN = 5
CXPB = 0.7
MUTPB = 0.2

def visualize_individuals(individuals, titles, layout):
    num_individuals = len(individuals)
    nrows, ncols = layout
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    
    if nrows * ncols < num_individuals:
        raise ValueError("Layout is too small for the number of individuals to display.")
    
    axes = np.array(axes).reshape(-1)  # Flatten in case of multiple rows
    
    for ax, individual, title in zip(axes, individuals, titles):
        antenna_array = np.array(individual[:NUM_ELEMENTS]).reshape(ARRAY_SHAPE)
        
        # Display the stored antenna elements directly without reprocessing
        ax.imshow(antenna_array, cmap='gray', interpolation='nearest')
        ax.set_title(f"{title}\n")
        ax.axis('off')
    
    # Hide any unused subplots
    for ax in axes[num_individuals:]:
        ax.axis('off')
    
    plt.tight_layout()
    plt.show()

def main():
    
    pop = toolbox.population(n=10)
    
    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    
    #pick two individuals from the population
    ind1 = pop[0]
    ind2 = pop[1]
    
    #copy both of them    
    copy_ind1 = toolbox.clone(ind1)
    copy_ind2 = toolbox.clone(ind2)
    
    #mate the two individuals
    child1, child2 = toolbox.mate(copy_ind1, copy_ind2)
    mutant1 = toolbox.mutate(toolbox.clone(child1))
    mutant2 = toolbox.mutate(toolbox.clone(child2))
    child12,child22 = toolbox.mate(copy_ind1, copy_ind2)    
    visualize_individuals([copy_ind1, child1, mutant1,child12, copy_ind2, child2, mutant2, child22],
                            ["Parent 1", "Child 1", "Mutant 1","Child 12", "Parent 2", "Child 2", "Mutant 2", "Child 22"],
                            layout=(2, 4))
    
    
# call main function
if __name__ == "__main__":
    random.seed(42)
    while True:
        main()
