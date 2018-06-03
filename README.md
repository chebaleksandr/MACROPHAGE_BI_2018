# Project "Investigation and development of the matabolic macrophage model"
## Students: А.А.Cheblockov, N.P.Rodina
## Supervisord: А.Н. Гайнуллина, А.А. Сергушичев

## Short project description

##### Macrophages - cells of the first line of immune defense: destroy pathogens (M1), maintain tissue homeostasis (M2).
##### The use of the metabolic FBA model allows one to see the coordination between the metabolic pathways at the level of the whole cell. However, the existing FBA model of macrophage metabolism has a number of inaccuracies, and therefore it does not reflect the latest ideas about their M1 activation, formulated in the course of molecular biological experiments.

### The aim of the project is to detect and fix the incorrectnesses of the FBA-model of macrophage metabolism.

## Short description of the used methods

##### Methods: 
* FBA (flux balance analysis): Under the conditions of our steady state, the result of multiplying the stoichiometric matrix, composed of all the cell responses, by the flow rate through these reactions is zero.
* FVA (Flux Variability Analysis): Method for evaluating the minimum and maximum range of each reaction flux that can still satisfy the constraints using a double LP problem (i.e. a maximization and a subsequent minimization) for each reaction of interest.

## Contents of the repository

##### This repository contains all the information relevant to the project's implementation as of March 31, 2013:
##### used models,
##### skripts,
##### description of the versions of the used soft

## Models used in the project 

##### Metabolic macrophage model (doi:10.1038/msb.2012.21)
##### RAW264_7_v2.xml
##### RAW264_7_v3.xml

## Skripts
#### IB_project_30_03_18.py (N.Rodina):
##### To run the script, you must download the required model, install all necessary libraries (see the Used programs section), and also change the path to the model for your computer in the script:

```
#reading the model
data_dir = cobra.test.data_dir
model = cobra.io.read_sbml_model(join(data_dir, "D:/Institute for bioinformatics/Project_sem2/RAW264_7_v2.xml"))
```
##### This script allows you to read the used model, as well as get all the information about it (the number of reactions, metabolites, etc.). Also, this script allows you to conduct FBA and FVA analysis. The results of the FVA analysis are output to the csv file.

```
Number of reactions in the model 1398
Number of metabolites in the model 1009
Number of genes in the model 768
<Solution 0.000 at 0xa1f91f3da0>
<Solution 0.062 at 0xa1f93b7e48>
```
##### Results of the FVA-analysis are written to a csv-file (out.csv).
#### fva_analysis.r (N.Rodina):
##### To run the script, you must download files (out.csv, MetaNetX_annot.txt, MetaNetX_analysis.txt, MetaNetX_gene_reaction.txt, met_lethal_sum1.txt), install all necessary libraries, and also set the path to the folder with all files.
##### This script allows you to read the the results of FVA analysis, annotate and visualize them.
#### Examples of visualization:
![All groups](https://github.com/chebaleksandr/MACROPHAGE_BI_2018/blob/master/all_rko.png)


#### Macrophage_Model.py (А. Cheblokov): 
##### To run the script, you must download the required model, install all necessary libraries (see the Used programs section), and also change the path to the model for your computer in the script:
```
mdl=cobra.io.read_sbml_model('/home/aleksandr/Downloads/RAW264_7_v3.xml')
```
##### In this script, the model RECON1.COMBINED.json is used to visualize the resulting solution

```
json_string = urllib.request.urlopen("https://raw.githubusercontent.com/escher/community-maps/master/RECON1/RECON1.COMBINED.json").read().decode('utf-8')
RECON = json.loads(json_string)
```

##### This script allows you to read the model from the xml file, and read the map from the json file, and also render with the Escher library. Also the function react_resemblance, allows you to compare the differences in the number of reactions between the card and the model.

## Used soft

##### Python 3.6.3 interpretator
##### Cobrapy library for Python: cobra (version 0.11.3) (https://cobrapy.readthedocs.io/en/latest/)
##### Pandas library for Python:  pandas (version 0.22.0)
##### Escher library for the visualization: Escher (version 1.6.0.)
##### A map for visualization with the use of Escher: RECON1.COMBINED.json (https://github.com/escher/community-maps/blob/master/RECON1/RECON1.COMBINED.json)
##### Numerical solutions soft (solvers): cplex, Gurobi Optimization

## Literature
##### 1. Aarash Bordbar et al., (2012)  Model-driven multi-omic data analysis elucidates metabolic immunomodulators of macrophage activation.  Molecular Systems Biology 8:558
##### 2. Natalie C. Duarte et al., (2006) Global reconstruction of the human metabolic network based on genomic and bibliomic data.  PNAS 104(6):1777–1782
##### 3. P. Kent Langston et al., (2017) Metabolism Supports Macrophage Activation. Front. Immunol. 8:61.
##### 4. Silvia Galván-Peña and Luke A. J. O’Neill (2014) Metabolic reprograming in macrophage polarization. Frontiers in Immunology, Inflammation 5:420:2
##### 5. Ryan DG, O’Neill LA. (doi: 10.1002/1873-3468.12744)
##### 6. Anna S. Blazier and Jason A. Papin. (2012) Integration of expression data in genome-scale metabolic network reconstructions. Front.      Physiol. 3:299
##### 7. Jeffrey D Orth et al., (2010) What is flux balance analysis? nature biotechnology 28(3): 245-248
##### 8. Gudmundsson, S., Thiele, I. Computationally efficient flux variability analysis. BMC Bioinformatics. 11, 489 (2010).





