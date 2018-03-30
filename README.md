# Проект "Изучение и разработка метаболической модели макрофагов"
## Студенты: А.А.Чеблоков, Н.П.Родина
## Руководители: А.Н. Гайнуллина, А.А. Сергушичев

## Краткое описание проекта

##### Макрофаги — клетки первой линии иммунной защиты: уничтожают патогены (М1), поддерживают тканевой гомеостазис (М2).
##### Использование метаболической FBA-модели позволяет увидеть координацию между метаболическими путями на уровне целой клетки. Однако существующая FBA-модель метаболизма макрофагов имеет ряд неточностей, в связи с чем она не отражает самые последние представления об их М1-активации, сформулированные в ходе молекулярно-биологических экспериментов.

### Целью настоящего проекта является обнаружить и исправить неточности FBA-модель метаболизма макрофагов.

## Краткое описание использованных методов

##### Метод: FBA (flux balance analysis): В условиях допущенного нами стационарного состояния результат умножения стехиометрической матрицы, составленной из всех реакций клетки, на скорость потоков через эти реакции равен нулю.

## Содержимое репозитория

##### Данный репозиторий содержит всю актуальную на 31.03.2018 информацию по выполнению проекта:
##### используемые модели,
##### скрипты,
##### описание версий используемых программ

## Использованные в работе модели 

##### Метаболическая модель макрофага (doi:10.1038/msb.2012.21)
##### RAW264_7_v2.xml
##### RAW264_7_v3.xml

## Скрипты
#### IB_project_30_03_18.py (Н.Родина):
##### Для запуска скрипта необходимо скачать требуемую модель, установить все необходимые библиотеки (см.раздел Использованные программы), а также изменить в скрипте путь к модели для своего компьютера:

```
#reading the model
data_dir = cobra.test.data_dir
model = cobra.io.read_sbml_model(join(data_dir, "D:/Institute for bioinformatics/Project_sem2/RAW264_7_v2.xml"))
```
##### Данный скрипт позволяет прочитать используемую модель, а также получить всю необходимую о ней информацию (количество реакций, метаболитов и др.). Также данный скрипт позволяет провести FBA и FVA анализ. Результаты FVA анализа выводятся в csv файл. 

```
Number of reactions in the model 1398
Number of metabolites in the model 1009
Number of genes in the model 768
<Solution 0.000 at 0xa1f91f3da0>
<Solution 0.062 at 0xa1f93b7e48>
```
##### Результаты FVA анализа выводятся в csv файл (out.csv).

##### Macrophage_Model.py (А. Чеблоков): 
##### Для запуска скрипта необходимо скачать требуемую модель, установить все необходимые библиотеки (см.раздел Использованные программы), а также изменить в скрипте путь к модели для своего компьютера:

```
mdl=cobra.io.read_sbml_model('/home/aleksandr/Downloads/RAW264_7_v3.xml')
```
##### В данном скрипте для визуализации полученного решения используется модель RECON1.COMBINED.json

```
json_string = urllib.request.urlopen("https://raw.githubusercontent.com/escher/community-maps/master/RECON1/RECON1.COMBINED.json").read().decode('utf-8')
RECON = json.loads(json_string)
```

##### Данный скрипт позволяет прочитать модель из xml файла, и прочитать карту из json файла, а также осуществить визуализацию при помощи библиотеки Escher. Также функция react_resemblance, позволяет сравнить различия в количестве реакций между картой и моделью.


## Использованные программы

##### Интерпретатор Python 3.6.3
##### Библиотека Cobrapy для Python: cobra (version 0.11.3) (https://cobrapy.readthedocs.io/en/latest/)
##### Библиотека Pandas для Python:  pandas (version 0.22.0)
##### Библиотека Escher для визуализации: Escher (version 1.6.0.)
##### Карта для визуализации с использованием Escher: RECON1.COMBINED.json (https://github.com/escher/community-maps/blob/master/RECON1/RECON1.COMBINED.json)
##### Программы для численного решения (солверы): cplex, Gurobi Optimization

## Литература
##### 1. Aarash Bordbar et al., (2012)  Model-driven multi-omic data analysis elucidates metabolic immunomodulators of macrophage activation.  Molecular Systems Biology 8:558
##### 2. Natalie C. Duarte et al., (2006) Global reconstruction of the human metabolic network based on genomic and bibliomic data.  PNAS 104(6):1777–1782
##### 3. P. Kent Langston et al., (2017) Metabolism Supports Macrophage Activation. Front. Immunol. 8:61.
##### 4. Silvia Galván-Peña and Luke A. J. O’Neill (2014) Metabolic reprograming in macrophage polarization. Frontiers in Immunology, Inflammation 5:420:2
##### 5. Ryan DG, O’Neill LA. (doi: 10.1002/1873-3468.12744)
##### 6. Anna S. Blazier and Jason A. Papin. (2012) Integration of expression data in genome-scale metabolic network reconstructions. Front.      Physiol. 3:299
##### 7. Jeffrey D Orth et al., (2010) What is flux balance analysis? nature biotechnology 28(3): 245-248







