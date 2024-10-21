from ifa_test import ifa_simulation
import random
from deap import base, creator, tools, algorithms

def evaluation_fun(individual):
    ifa_h, ifa_l, ifa_fp = individual

    Sim_CSX = 'IFA.xml'
    showCad = False
    post_proc_only = False
    unit = 1e-3
    substrate_width = 21
    substrate_length = 83.15
    substrate_thickness = 1.5
    gndplane_position = 0
    substrate_cells = 4
    #ifa_h = 5.586
    #ifa_l = 18.0
    ifa_w1 = 0.6
    ifa_w2 = 0.6
    ifa_wf = 0.6
    #ifa_fp = 1.108
    ifa_e = 0.5
    substrate_epsR = 4.5
    feed_R = 50
    min_freq = 2.4e9
    center_freq = 2.45e9
    max_freq = 2.5e9
    min_size = 0.2,  # minimum automesh size
    plot = False

    freq, s11_dB, f_res, Zin, P_in = ifa_simulation(Sim_CSX=Sim_CSX,
                                                   showCad=showCad,
                                                   post_proc_only=post_proc_only,
                                                   unit=post_proc_only,
                                                   substrate_width=substrate_width,
                                                   substrate_length=substrate_length,
                                                   substrate_thickness=substrate_thickness,
                                                   gndplane_position=gndplane_position,
                                                   substrate_cells=substrate_cells,
                                                   ifa_h=ifa_h,
                                                   ifa_l=ifa_l,
                                                   ifa_w1=ifa_w1,
                                                   ifa_w2=ifa_w2,
                                                   ifa_wf=ifa_wf,
                                                   ifa_fp=ifa_fp,
                                                   ifa_e=ifa_e,
                                                   substrate_epsR=substrate_epsR,
                                                   feed_R=feed_R,
                                                   min_freq=min_freq,
                                                   center_freq=center_freq,
                                                   max_freq=max_freq,
                                                   plot=plot,
                                                   min_size=min_size)



if __name__ == "__main__":

    # Genetic Algorithm setup
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    toolbox.register("ifa_h", random.uniform, 10, 50)
    toolbox.register("ifa_l", random.uniform, 2, 10)
    toolbox.register("ifa_fp", random.uniform, 1, 10)
    toolbox.register("individual", tools.initCycle, creator.Individual,
                     (toolbox.attr_length, toolbox.attr_width, toolbox.attr_feed), n=1)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", evaluation_fun)
    toolbox.register("mate", tools.cxBlend, alpha=0.5)
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=0.2)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # Genetic Algorithm execution
    population = toolbox.population(n=30)
    NGEN = 40
    for gen in range(NGEN):
        offspring = algorithms.varAnd(population, toolbox, cxpb=0.5, mutpb=0.2)
        fits = map(toolbox.evaluate, offspring)
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit
        population = toolbox.select(offspring, k=len(population))
    best_ind = tools.selBest(population, 1)[0]

    print("Best individual is %s with fitness %s" % (best_ind, best_ind.fitness.values))
