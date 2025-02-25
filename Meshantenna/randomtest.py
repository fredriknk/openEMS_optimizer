import numpy as np
from scipy.ndimage import gaussian_filter, zoom
from deap import base, creator, tools
import matplotlib.pyplot as plt

# -----------------------------
# Global parameters
# -----------------------------
ARRAY_SIZE = (10, 20)

# -----------------------------
# 1) Define the clustered array generator
# -----------------------------
def generate_clustered_array(array_size=(40, 40), mean=0.5, std_dev=0.2, correlation_length=5, threshold=0.509):
    random_field = np.random.normal(mean, std_dev, array_size)
    smoothed_field = gaussian_filter(random_field, sigma=correlation_length)
    antenna_array = smoothed_field > threshold
    antenna_array = antenna_array.astype(int)
    return antenna_array.flatten()

def generate_random_array(array_size=(40, 40), threshold=0.4,must_be_connected=5):
    connected = False
    while not connected:
        random_field = np.random.rand(*array_size)
        antenna_array = random_field > threshold
        antenna_array = antenna_array.astype(int)
        if antenna_array[-1,must_be_connected]:
            connected = True
    return antenna_array.flatten()

# -----------------------------
# 2) Define the custom individual initialization
# -----------------------------
def init_clustered_individual():
    flat_antenna_array = generate_random_array(array_size=ARRAY_SIZE)
    return creator.Individual(flat_antenna_array.tolist())

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("individual", init_clustered_individual)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# Use our custom mate operator
toolbox.register("mate", tools.cxUniform, indpb=0.5)

# Use a small mutation
toolbox.register("mutate", tools.mutFlipBit, indpb=0.01)

# Basic selection
toolbox.register("select", tools.selTournament, tournsize=3)

# A simple evaluation: count the number of 1s
def evaluate(ind):
    return (sum(ind),)
toolbox.register("evaluate", evaluate)

# -----------------------------
# 6) Demonstrate creation, mating, mutation, and visualization
# -----------------------------
pop = toolbox.population(n=2)  # 2 parents

# Evaluate
for ind in pop:
    ind.fitness.values = toolbox.evaluate(ind)

parent1, parent2 = pop[0], pop[1]

# Clone them as "offspring"
offspring1 = toolbox.clone(parent1)
offspring2 = toolbox.clone(parent2)

# Apply our custom mate (crossover) in-place
toolbox.mate(offspring1, offspring2)
del offspring1.fitness.values, offspring2.fitness.values

# Apply a little mutation
toolbox.mutate(offspring1)
toolbox.mutate(offspring2)

# Evaluate offspring
for child in (offspring1, offspring2):
    child.fitness.values = toolbox.evaluate(child)

# Print fitness
print("Parent1 Fitness:", parent1.fitness.values)
print("Parent2 Fitness:", parent2.fitness.values)
print("Offspring1 Fitness:", offspring1.fitness.values)
print("Offspring2 Fitness:", offspring2.fitness.values)

# Visualize
fig, axs = plt.subplots(2, 2, figsize=(8,8))

axs[0,0].imshow(np.array(parent1).reshape(ARRAY_SIZE), cmap='Greys')
axs[0,0].set_title(f"Parent1: {parent1.fitness.values[0]}")

axs[0,1].imshow(np.array(parent2).reshape(ARRAY_SIZE), cmap='Greys')
axs[0,1].set_title(f"Parent2: {parent2.fitness.values[0]}")

axs[1,0].imshow(np.array(offspring1).reshape(ARRAY_SIZE), cmap='Greys')
axs[1,0].set_title(f"Offspring1: {offspring1.fitness.values[0]}")

axs[1,1].imshow(np.array(offspring2).reshape(ARRAY_SIZE), cmap='Greys')
axs[1,1].set_title(f"Offspring2: {offspring2.fitness.values[0]}")

print(np.array(offspring2).reshape(ARRAY_SIZE))

plt.tight_layout()
plt.show()
