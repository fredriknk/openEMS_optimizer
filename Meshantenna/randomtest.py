import numpy as np
from scipy.ndimage import gaussian_filter
from deap import base, creator, tools
import matplotlib.pyplot as plt

# Step 1: Define the clustered array generation function
def generate_clustered_array(array_size=(40, 40), mean=0.5, std_dev=0.3, correlation_length=2, threshold=0.501):
    random_field = np.random.normal(mean, std_dev, array_size)
    smoothed_field = gaussian_filter(random_field, sigma=correlation_length)
    antenna_array = smoothed_field > threshold
    antenna_array = antenna_array.astype(int)
    flat_antenna_array = antenna_array.flatten()
    return flat_antenna_array

# Step 2: Define the custom individual initialization function
def init_clustered_individual():
    flat_antenna_array = generate_clustered_array()
    return creator.Individual(flat_antenna_array.tolist())

# Step 3: Setup DEAP toolbox
creator.create("FitnessMax", base.Fitness, weights=(1.0,))  # Adjust weights as needed
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("individual", init_clustered_individual)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# Register genetic operators
toolbox.register("mate", tools.cxOnePoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)

# Step 4: Define the evaluation function
def evaluate(individual):
    # Implement your evaluation logic here
    # For example, count the number of antenna elements
    num_antennas = sum(individual)
    # Return a tuple as required by DEAP
    return (num_antennas,)  # Maximize the number of antennas as an example

toolbox.register("evaluate", evaluate)

# Step 5: Test the initialization and population generation
population_size = 100
population = toolbox.population(n=population_size)

# Evaluate the initial population
fitnesses = list(map(toolbox.evaluate, population))
for ind, fit in zip(population, fitnesses):
    ind.fitness.values = fit

# Print and visualize the first individual
print("First individual:")
print(population[0])
print("Fitness:", population[0].fitness.values)

# Visualize the first individual's antenna array
antenna_array_2d = np.array(population[0]).reshape((40, 40))

plt.imshow(antenna_array_2d, cmap='Greys', interpolation='none')
plt.title('First Individual - Clustered Antenna Array')
plt.show()
