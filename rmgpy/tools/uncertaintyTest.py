################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import unittest
import os
from rmgpy.tools.uncertainty import *
import rmgpy
from rmgpy.data.rmg import RMGDatabase 
        
class uncertaintyTest(unittest.TestCase):

    def testUncertainties(self):
        """
        Test that the parameter uncertainty engines can be loaded
        """
        rmgpyName = os.path.dirname(rmgpy.__file__)
        dataDir = os.path.join(rmgpyName,'test_data','testing_database')
        chemDir = os.path.join(rmgpyName,'test_data','parsing_data')
        chemkinFile = os.path.join(chemDir,'chem_annotated.inp')
        spcDict = os.path.join(chemDir,'species_dictionary.txt')
        
        uncertainty = Uncertainty(outputDirectory='chemDir')
        uncertainty.loadModel(chemkinFile, spcDict)
        
        #load database properly
        
        kineticsFamilies = 'all'
        
        kineticsDepositories = ['training']

        thermoLibraries = ['primaryThermoLibrary']

        reactionLibraries = ['GRI-Mech3.0']
            
        uncertainty.database = RMGDatabase()
        uncertainty.database.load(dataDir, 
                          kineticsFamilies=kineticsFamilies, 
                          kineticsDepositories=kineticsDepositories,
                          thermoLibraries=thermoLibraries,
                          reactionLibraries=reactionLibraries,
                          )
        
        # Prepare the database by loading training reactions but not averaging the rate rules
        for familyLabel, family in uncertainty.database.kinetics.families.iteritems():
            family.addKineticsRulesFromTrainingSet(thermoDatabase=uncertainty.database.thermo)
        
            family.fillKineticsRulesByAveragingUp(verbose=True)
            
        uncertainty.extractSourcesFromModel()
        uncertainty.compileAllSources()
        uncertainty.assignParameterUncertainties()
        gParamEngine = ThermoParameterUncertainty()
        kParamEngine = KineticParameterUncertainty()
        