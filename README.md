################################################################################
################################################################################

This repository contains data from:

################################################################################

Carrion decomposition in a subtropical biodiversity experiment

dx.doi.org/10.6084/m9.figshare.30444791

Finn Rehling1, Matteo Dadda2, Marc Nagel1, Georg Albert3, Helge Bruelheide4,5, Jing-Ting Chen3,6, Heike Feldhaar2, Felix Fornoff1, Arong Luo6, 
Massimo Martini1,6, Xiao-Yu Shi6, Michael Staab7, Xianglu Deng8, Xiaojuan Liu8, Qing-Song Zhou6, Chao-Dong Zhu6, Alexandra-Maria Klein1.

1Chair of Nature Conservation and Landscape Ecology, Albert-Ludwigs-University Freiburg, Freiburg, Germany.

2Animal Population Ecology, Bayreuth Center of Ecology and Environmental Research (BayCEER), University of Bayreuth, Bayreuth, Germany.

3Forest Nature Conservation, Georg-August-University Göttingen, Göttingen, Germany.

4Institute of Biology/Geobotany and Botanical Garden, Martin Luther University Halle-Wittenberg, Halle (Saale), Germany.

5German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena-Leipzig, Leipzig, Germany.

6State Key Laboratory of Animal Biodiversity Conservation and Integrated Pest Management, Institute of Zoology, Chinese Academy of Sciences, Beijing, China.

7Animal Ecology and Trophic Interactions, Institute of Ecology, Leuphana University Lüneburg, Lüneburg, Germany.

8Institute of Botany, Chinese Academy of Sciences, Beijing, China.

#################################################################################

Data holder: Finn Rehling
Chair of Nature Conservation and Landscape Ecology, Albert-Ludwigs-University Freiburg, Freiburg, Germany.

################################################################################


The package contains two data file (in .csv format) that were used in the analyses of the above-mentioned publication.

The structures of the data files are listed below with detailed descriptions of their contents.


################################################################################

1. Decomposer diversity.csv
Number of carrion-associated arthropods in carrion-baited traps in BEF-China (site A and B); sampled in May/June 2023 and July 2024.

 $ site             : fac, experimental site (A or B)
 
 $ plot          : fac, experimental plot nested in site 
 
 $ year      : num, study year (2023 or 2024)
 
 $ fly.1          : int, abundance of morphospecies 'fly.1'
 
 $ fly.3:formicidae.k         : int, abundance of other morphospecies
 
 $ vespidae.1:mantidae.y6         : int, abundance of predatory morphospecies; not included in the analysis of decomposers


################################################################################

> Carrion decomposition.csv

Data on carrion decomposition and associated arthropods in BEF-China (site A and B); sampled in May/June 2023 and July 2024.

$ site : fac, experimental site (A or B)

$ plot : fac, experimental plot nested in site

$ canopy.cover : num, canopy cover of the forest at plot level

$ slope.steepness : num, slope steepness of the forest at plot level

$ tree.richness : num, tree richness of the forest at plot level

$ location : fac, location of carrion nested in plot

$ x : fac, x coordinate of carrion location nested in plot

$ y : fac, y coordinate of carrion location nested in plot

$ tree.species : fac, tree species around which the carrion was placed

$ year : fac, study year (2023 or 2024)

$ ID : fac, ID for carrion nested in location/plot/site

$ initial.date : num, date when the experiment started (dd-mm-yyyy)

$ date.photo : date, date when the carrion photo was taken (dd-mm-yyyy)

$ days : num, number of days until the carrion photo was taken (duration of decomposition)

$ days.fac : fac, "early" = first photo, "late" = second photo

$ initial.mass : num, initial fresh mass of carrion [g * 10]

$ remain.mass : num, remaining fresh mass of carrion after two days [g * 10]

$ removal : fac, carrion removed by scavengers (0 = no, 1 = yes)

$ deco.raw : ord, averaged decomposition score (1 = "fresh", 7 = "remains")

$ deco.adj : ord, adjusted decomposition score for differences in observation date (1–7)

$ deco.EN : ord, decomposition score by observer EN (1–7)

$ deco.MN : ord, decomposition score by observer MN (1–7)

$ deco.FR : ord, decomposition score by observer FR (1–7)

$ deco.JC : ord, decomposition score by observer JC (1–7)

$ rel.remain : num, fresh remaining mass of carrion [%] after two days

$ remain : num, predicted absolute remaining mass of carrion [g * 10] after two days

$ loss : num, predicted absolute mass lost [g * 10] after two days

$ prime.deco : fac, primary decomposer ("ants", "flies", "both", "unknown")

$ fly.prob : fac, flies as primary decomposer (0 = no, 1 = yes)

$ ant.prob : fac, ants as primary decomposer (0 = no, 1 = yes)

$ both.prob : fac, ants and flies as primary decomposer (0 = no, 1 = yes)

$ burried : fac, burried by ants (0 = no, 1 = yes, 2 = half-burried)

$ no.flies : num, number of carrion-decomposing flies in trap at plot level

$ flies_Estimator : num, species richness of flies in trap at plot level

$ no.other : num, number of other carrion-decomposing arthropods in trap at plot level

$ other_Estimator : num, species richness of other arthropods in trap at plot level

$ no.ants : num, number of carrion-decomposing ants in trap at plot level

$ ants_Estimator : num, species richness of ants in trap at plot level
