# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:11:10 2017

@author: AB Sanyal
"""

import numpy as np
import random
import matplotlib.pyplot as plt
import time

#Initial values
dna_len = 50
desired_x, desired_y = 20, 0
family_size = 20
mutation_rate = 0.05

fitness_list = []
gen_list = []

#Dictionary of amino-acids
bases = dict([("0", "U"), ("1", "D"), ("2", "L"), ("3", "R")])

#Show a formatted sequence where elements are separated by a given separator
def show_seq(seq, separator):
    i = 0
    p = ""
    while (i < len(seq)):
        p += str(seq[i])
        if (i < len(seq) - 1):
            p += str(separator)
        i += 1
    return p

#Generate a random DNA of a particular length n
def generate_dna(n):
    i = 0
    p = []
    while (i < n):
        r = random.randint(0, 3)
        p.append(r)
        i += 1
    return p

#Norm of a coordinate
def norm(coordinate):
    p = 0
    for component in coordinate:
        p += component ** 2
    return p ** 0.5

#Create an easy to read form
def readable_form(dna):
    p = []
    for base in dna:
        p.append(bases[str(base)])
    return p

#Read the DNA and give the final deviation from desired location
def find_final_deviation(dna):
    x, y = 0, 0
    for base in dna:
        if base == 0: y += 1
        if base == 1: y += -1
        if base == 2: x += -1
        if base == 3: x += 1
    return norm((x - desired_x, y - desired_y))

#Point Mutation
def mutate(dna):
    mutated_dna = []
    for base in dna:
        p = random.uniform(0, 1)
        if p < mutation_rate:
            basenew = random.choice([0, 1, 2, 3])
            while (basenew == base):
                basenew = random.choice([0, 1, 2, 3])
        else:
            basenew = base
        mutated_dna.append(basenew)
    return mutated_dna


#Cross-breed at random cross-over point
def crossover(p1, p2):
    cop = np.random.randint(1, dna_len)
    child = p1[:cop] + p2[cop:]
    child = mutate(child)
    return child

#Sexually recombined offsprings
def create_recombinants(p1, p2, n):
    i = 1
    family = []
    while (i <= n):
        if (np.random.uniform(0, 1) <= 0.5):
            child = crossover(p1, p2)
        else:
            child = crossover(p2, p1)
        family.append(child)
        i += 1
    return family

#Create a family of n random paths
def create_family(n):
    family = []
    i = 0
    while (i < n):
        family.append(generate_dna(dna_len))
        i += 1
    return family

#Main code
family = []

current_best_fitness = 10**6

best_gene = []
random_gene = []

#Create a pioneer family
family = create_family(family_size)

gen_no = 0 #Gen 0 is the F0

while (current_best_fitness != 0):
    #Check for the fittest
    for this_gene in family:
        this_gene_fitness = find_final_deviation(this_gene)
    if this_gene_fitness <= current_best_fitness:
        current_best_fitness = this_gene_fitness
        best_gene = this_gene

    #Print some stuff every 1000 generations
    if (gen_no < 1000 or gen_no % 1000 == 0):
        print("Generation No:", str(gen_no), "Best gene's fitness:", str(1/(1 + current_best_fitness)))
        print("The best gene is:")
        print(show_seq(readable_form(best_gene),""))
        print(":" * 100)

    #Plot path of best offspring
    xcor_list = []
    ycor_list = []
    x = 0
    y = 0
    for base in best_gene:
        if base == 0: y += 1
        if base == 1: y += -1
        if base == 2: x += -1
        if base == 3: x += 1
        xcor_list.append(x)
        ycor_list.append(y)
    plt.plot(xcor_list, ycor_list)

    #fitness of best offspring
    fitness_list.append(1/(1 + current_best_fitness))
    gen_list.append(gen_no)

    #Select the best and a random parent
    p1 = best_gene
    p2 = random.choice(family)
    while p1 == p2: #Prevent self-cross
        p2 = random.choice(family)

    #Create offsprings
    family = create_recombinants(p1, p2, family_size)

    #Advance the generation counter
    gen_no += 1

    #time.sleep(0.5)

#Print info about final gen
print("FINAL GENERATION")
print("Generation No:", str(gen_no - 1), "Best gene's fitness:", str(1/(1 + current_best_fitness)))
print("The best gene is:")
print(show_seq(readable_form(best_gene),""))
print(":" * 100)

#Plot the best path
xcor_list = []
ycor_list = []
x = 0
y = 0
for base in best_gene:
    if base == 0: y += 1
    if base == 1: y += -1
    if base == 2: x += -1
    if base == 3: x += 1
    xcor_list.append(x)
    ycor_list.append(y)

plt.plot(xcor_list, ycor_list, '-o')
plt.show()

plt.plot(gen_list, fitness_list)

plt.show()