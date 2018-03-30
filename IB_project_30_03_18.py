#Project 2 semester
#created by Natalia Rodina
#last update 30.03.2018

#This skirpt opens the macrophage model
import cobra
import cobra.test
import os
import lxml
from os.path import join
from cobra import Reaction, Metabolite, Model
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.flux_analysis import pfba
import pandas
import matplotlib.pyplot as plt
#import plot
from cobra.flux_analysis import flux_variability_analysis

#import plot_helper


#reading the model
data_dir = cobra.test.data_dir
model = cobra.io.read_sbml_model(join(data_dir, "D:/Institute for bioinformatics/Project_sem2/RAW264_7_v2.xml"))

#shows the number of reactions, metabolites and genes in the model
print('Number of reactions in the model', len(model.reactions))
print('Number of metabolites in the model', len(model.metabolites))
print('Number of genes in the model', len(model.genes))
#print('Model', model)

#working with one reaction
#print ('Reaction number 2', model.reactions[29])
#react_1 = model.reactions.get_by_id("10FTHF5GLUtl")
#print(react_1.name)
#print(react_1.reaction)
#print(react_1.lower_bound, "< pgi <", react_29.upper_bound)
#print(react_1.reversibility)


#looples solution
nominal = model.optimize()
loopless = loopless_solution(model)
print(loopless)
df = pandas.DataFrame(dict(loopless=loopless. fluxes, nominal=nominal.fluxes)) #as dataframe
df.plot.scatter(x='loopless', y='nominal') #plot the solution
plt.show()

solution = model.optimize() #optimize model
print(solution)
#model.solver = 'glpk'

#FVA
print(flux_variability_analysis(model, model.reactions, fraction_of_optimum=0.9))
dataframe = flux_variability_analysis(model, model.reactions) #results to a dataframe
dataframe.to_csv('out.csv') #writes the dataframe to csv file
#print (cobra.flux_analysis.flux_variability_analysis(
#    model, model.reactions[:10], fraction_of_optimum=0.9))


#model.optimize()
#model.summary(fva=0.95)
