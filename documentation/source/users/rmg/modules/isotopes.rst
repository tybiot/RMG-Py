.. _isotopes:

********
Isotopes
********

Isotopic enrichment can be indicated in a molecular structure's adjacency list. 
The example below is methane with an isotopically labeled carbon of isotope number 
13, which is indicated with ``i13``::


    1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 H u0 p0 c0 {1,S}
    5 H u0 p0 c0 {1,S}

The less resource intensive isotope generation algorithm involves a two step process. 
A model without any isotopes is generated, and then it is used to create a model with 
isotopes by labeling all molecules and then generating reactions between species. 
This requires two input files. The first one is a standard RMG input file without 
isotopes. The second input file can contain an optional parameter `maximumIsotopicAtoms` 
in the `options` block. This limits how many isotopic labels can exist in one 
molecule, and can reduce the model size complexity. 

The algorithm can be run with the command::

    python \$rmg/scripts/isotopes.py path/to/input.py

The input file in this case should be free of isotopes, since the model will 
generate them automatically. One extra species contraint is available for 
isotope runs, \lstinline{maximumIsotopicAtoms}, which restricts the number 
of isotopes of a particular element that can be labeled on a molecule. For 
low-enrichement scenarios, this could reduce computational costs.

The script generates a mechanism without isotopes first, finds all the 
isotopomers that fit within the ``maximumIsotopicAtoms`` constraint, and then 
generates the reactions between all labeled molecules. If you already have a 
model (which includes atom mapping in RMG's format) which you would like to 
add isotope labels to, you can use the command::

    python \$rmg/scripts/isotopes.py path/to/input.py  --original path/to/model/directory

The desired model input should be stored within with structure 
``chemkin/chem_annotated.inp` and ``chemkin/species_dictionary.txt``.

Some functionality in RMG software conflicts with how models are generated. 
Reaction libraries need to explicitly mention the reactions of isotopomers 
so they match various isotopes. If kinetic isotope effects for a specific 
reaction are to be included, they could be input in a reaction library, but 
this would need to be done manually. 

Following the generation, a number of diagnostics are ran to check model 
accuracy. Isotopomers are checked to ensure their symmetries are consistent. 
Then, the reaction path degeneracy among reactions differing only in isotope 
labeling is checked to ensure it is consistent with the symmetry values of reactions. 
If one of these checks throws a warning, the model will likely excibit non-natural 
fluctuations in enrichment.

This algorithm is currently limited to systems without significant kinetic 
isotope effects, and currently only works for Carbon-13.