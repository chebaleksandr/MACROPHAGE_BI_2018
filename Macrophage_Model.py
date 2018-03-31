#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
autor: Cheblokov_AA
"""
import cobra
import cobra.test
import urllib
from os.path import join
import libsbml
import urllib
import json
import escher


# Записываем карту реакций
json_string = urllib.request.urlopen("https://raw.githubusercontent.com/escher/community-maps/master/RECON1/RECON1.COMBINED.json").read().decode('utf-8')
RECON = json.loads(json_string)


# Записываем модель
model=cobra.io.read_sbml_model('/home/aleksandr/Downloads/RAW264_7_v3.xml')

#Обращение к реакциям, метаболитам, генам внутри модели
print(len(model.reactions)) 
print(len(model.metabolites))
print(len(model.genes))

model.reactions.PGI.bounds # номер реакции, затем обращение к переменной с ограничениями

pgi = model.reactions.get_by_id("PGI") # Записываем в переменную реакцию

print(pgi.name) # имена и реакция
print(pgi.reaction)
print(pgi.lower_bound, "< pgi <", pgi.upper_bound)
print(pgi.reversibility)

pgi.check_mass_balance()

pgi.add_metabolites({model.metabolites.get_by_id("h_c"): -1}) # Добавление метаболита
pgi.reaction

pgi.subtract_metabolites({model.metabolites.get_by_id("h_c"): -1})
# pgi.subtract_metabolites \ add_metabolites({Что: С каким коэффициентом (место и количество)})


# Функция для сравнения модели и карты

def react_resemblance(md_1, map_2):
    # Variables
    l_1 = md_1.reactions
    l_2 = []
    for key in map_2[1]['reactions'].keys():
        l_2.append(map_2[1]['reactions'][key]['bigg_id'])
    
    list_1 = []
    list_2 = []
    iter_list_1=[]
    iter_list_2=[]
    global pos_list_1, pos_list_2
    pos_list_1, pos_list_2=[],[]

    i=0
    for element in l_1:
        iter_list_1.append(i)
        pos_list_1.append(i)
        list_1.append(element.id)
        i+=1
    i=0
    for element in l_2:
        iter_list_2.append(i)
        pos_list_2.append(i)
        list_2.append(element)
        i+=1
        
    for each in iter_list_1:
        for every in iter_list_2:
            if list_1[each]==list_2[every]:
                print(each, every)
                if each in pos_list_1:
                    pos_list_1.remove(each)
                if every in pos_list_2:
                    pos_list_2.remove(every)
    print('Model, availability of reactions №:',pos_list_1,'\n'+'='*40+'\nMap, availability of reactions №:',pos_list_2)
    

react_resemblance(model,RECON)

# Список реакций на карте
list_of_reactions = []
for key in RECON[1]['reactions'].keys():
    list_of_reactions.append(RECON[1]['reactions'][key]['bigg_id'])
    

# Вывод реакций различных для карты и модели
with open('/home/aleksandr/Documents/result_1.txt','w') as txt:
    txt.write('-'*50+'\n\tEXCESS REACTIONS IN MODEL\n'+'-'*50+'\n')
    for i in pos_list_1:
        txt.write(str(model.reactions[i].id)+"\n")   
    txt.write('-'*50+'\n\tEXCESS REACTIONS IN RECON MAP\n'+'-'*50+'\n')
    for i in pos_list_2:
        txt.write(list_of_reactions[i]+"\n")


# Удаление из карты реакций, не соотносящихся с моделью
list_of_keys=[]
for key in RECON[1]['reactions'].keys():
    list_of_keys.append(key)
    
delete_list=[]
for i in pos_list_2:
    delete_list.append(list_of_keys[i])

for deletion in delete_list:
    RECON[1]['reactions'].pop(deletion)

len(RECON[1]['reactions'])

with open("/home/aleksandr/Documents/NEW_1.json",'w') as out:
    json.dump(RECON,out)


# Создание карты с моделью
model_map = escher.Builder(map_json = 'https://raw.githubusercontent.com/escher/community-maps/master/RECON1/RECON1.COMBINED.json',\
                      model = model)
model_map_2 = escher.Builder(map_json = "/home/aleksandr/Documents/NEW_1.json",\
                      model = model)

# Отрисовка карты в браузере
model_map.display_in_browser()
model_map_2.display_in_browser()
