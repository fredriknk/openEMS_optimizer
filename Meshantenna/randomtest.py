import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from deap import base, creator, tools, algorithms

# Define the shape of the 2D binary array
ARRAY_SHAPE = (40, 40)

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


def generate_clustered_array(array_size=(40, 40), std_dev_range = (0.2, 0.25),  correlation_length_range = (2, 5.0), threshold_range = (0.49, 0.52)):
    """
    Generate a clustered binary array using Gaussian smoothing.
    """
    mean = 0.5
    std_dev = random.uniform(*std_dev_range)
    correlation_length = random.uniform(*correlation_length_range)
    threshold = random.uniform(*threshold_range)
    
    random_field = np.random.normal(mean, std_dev, array_size)
    smoothed_field = gaussian_filter(random_field, sigma=correlation_length)
    antenna_array = smoothed_field > threshold
    antenna_array = antenna_array.astype(int)
    flat_antenna_array = antenna_array.flatten()
    return flat_antenna_array.tolist()

def init_individual():
    """
    Initialize an individual with a clustered antenna array.
    """
    antenna_elements = generate_clustered_array(
        array_size=ARRAY_SHAPE
    )
    return creator.Individual(antenna_elements)

toolbox.register("individual", init_individual)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

from scipy.ndimage import gaussian_filter
import numpy as np

def smoothed_mean_crossover(parent1, parent2, sigma=1.0, threshold=0.5):
    """
    Perform smoothed mean crossover between two parents.
    Applies Gaussian blur to both, averages them, thresholds to produce binary children,
    and returns two new child individuals without modifying the parents.
    
    Args:
        parent1 (list): First parent individual (list of 0s and 1s).
        parent2 (list): Second parent individual (list of 0s and 1s).
        sigma (float): Standard deviation for Gaussian blur.
        threshold (float): Threshold for binarization.
        
    Returns:
        tuple: Two new child individuals (lists of 0s and 1s).
    """
    # Ensure that ARRAY_SHAPE is defined globally or pass it as a parameter
    global ARRAY_SHAPE  # Replace with parameter if preferred
    
    # Reshape parents into 2D arrays
    array1 = np.array(parent1).reshape(ARRAY_SHAPE)
    array2 = np.array(parent2).reshape(ARRAY_SHAPE)
    
    # Apply Gaussian blur to both arrays
    r = random.random()
    ri = 1.-r
    
    blurred1 = gaussian_filter(array1.astype(float), sigma=sigma)*r
    blurred2 = gaussian_filter(array2.astype(float), sigma=sigma)*ri
    
    r1 = random.random()
    ri1 = 1.-r1
    
    blurred3 = gaussian_filter(array1.astype(float), sigma=sigma)*r1
    blurred4 = gaussian_filter(array2.astype(float), sigma=sigma)*ri1
    
    # Average the blurred arrays
    averaged = (blurred1 + blurred2) / 2.0
    averaged2 = (blurred3 + blurred4) / 2.0
    # Threshold the averaged array to obtain binary children
    binary_child = (averaged > threshold*random.random()).astype(int)
    binary_child2 = (averaged2 > threshold*random.random()).astype(int)
    # Flatten the 2D array back to a 1D list
    child1 = binary_child.flatten().tolist()
    child2 = binary_child2.flatten().tolist()  # Both children are identical in this approach
    
    # Alternatively, to introduce slight variations between children, you can add small noise
    # Uncomment the following lines to apply random perturbations to child2
    # noise = np.random.normal(0, 0.05, binary_child.shape)
    # averaged_with_noise = averaged + noise
    # binary_child2 = (averaged_with_noise > threshold).astype(int)
    # child2 = binary_child2.flatten().tolist()
    
    return child1, child2

def block_mutation(individual, block_size=(5,5), indpb=0.1):
    """
    Perform block-based mutation on an individual.
    """
    # Determine the number of blocks
    n_blocks_row = ARRAY_SHAPE[0] // block_size[0]
    n_blocks_col = ARRAY_SHAPE[1] // block_size[1]
    
    # Choose a random block to mutate
    block_row = random.randint(0, n_blocks_row - 1)
    block_col = random.randint(0, n_blocks_col - 1)
    
    # Calculate start and end indices
    start_idx = block_row * block_size[0] * ARRAY_SHAPE[1] + block_col * block_size[1]
    end_idx = start_idx + block_size[0] * ARRAY_SHAPE[1]
    
    # Flip bits in the block with probability indpb
    for i in range(start_idx, end_idx):
        if random.random() < indpb:
            individual[i] = 1 - individual[i]
    
    return (individual,)

toolbox.register("mate", smoothed_mean_crossover)
toolbox.register("mutate", block_mutation, block_size=BLOCK_SIZE, indpb=0.2)
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

    visualize_individuals([ind1, ind2, child1, child2],
                            ["Parent 1", "Parent 2", "Child 1", "Child 2"],
                            layout=(2, 2))
    
# call main function
if __name__ == "__main__":
    random.seed(42)
    while True:
        main()
