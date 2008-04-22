//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include "G4CompositeEMDataSet.hh"
#include "G4ShellEMPDataSet.hh"
#include "G4EMPDataSet.hh"
#include "G4VEMPDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"

#include "AIDA/AIDA.h"

int main()
{
  
    AIDA::ITree* tree;
    AIDA::IAnalysisFactory* analysisFactory;
    AIDA::ITreeFactory* treeFactory;

    AIDA::IPlotterFactory* plotterFactory;
    AIDA::IPlotter* plotter;

    AIDA::IDataPointSet* energ;
    AIDA::IDataPointSet* cross;
    AIDA::IDataPointSetFactory* curvFactory;

    analysisFactory = AIDA_createAnalysisFactory();
    treeFactory = analysisFactory->createTreeFactory();
    tree = treeFactory->create("CrossSection.xml");
    curvFactory = analysisFactory->createDataPointSetFactory(*tree);
 
    std::vector<G4double> engVal; 
    std::vector<G4double> crosVal;         
    std::vector<G4double> engErr; 
    std::vector<G4double> crosErr; 

  G4int z1;
  G4int z2; 
  G4double energy1;
  G4String fileName; 

 G4cout <<"----------------------------------------------------------------------------" << G4endl; 
 G4cout <<"Enter atomic number of the incident particles (protons or alpha): " << G4endl;
 G4cout <<"----------------------------------------------------------------------------" << G4endl;
 G4cin >> z1;

 if (z1 == 1)
   { fileName = "kcsPaul/kcs-";} 
  else
    {
      if (z1 == 2)
	{ fileName = "kacsPaul/kacs-";}
	
    }

  
 std::vector<G4double> energies;

   for (G4double i=-0.5; i<4.; i=i+0.1)
     {
       energies.push_back(std::exp(i));
     } 


 G4cout <<""<<G4endl;
 G4cout <<"----------------------------------------------------------------------------" << G4endl; 
 G4cout <<"Enter atomic number of element you want to check its K_Shell CrossSection: " << G4endl;
 G4cout <<"----------------------------------------------------------------------------" << G4endl;
 G4cin >> z2;

  
  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();

  G4VEMPDataSet* dataSet;

  dataSet = new G4EMPDataSet(z2,interpolation);
  
  dataSet->LoadData(fileName);
   
 
 for (size_t k=0; k<energies.size();k++)                        //Cross section for each incident energy
	{

	  energy1 = energies[k]/MeV;

          G4double sigma = dataSet->FindValue(energy1) / barn;
      
	 G4cout <<"with incident particles of "<<energy1/MeV <<" MeV, the K_CrossSection, of the element Z=" <<z2 <<", is "<<sigma <<" barn" << G4endl;
	 G4cout <<"-----------------------------------------------------------------------------------------" << G4endl; 
         
          engVal.push_back(energy1/MeV);
          crosVal.push_back(sigma);
          engErr.push_back(0);
          crosErr.push_back(0);
	  
	   }  

  delete dataSet;

      
       plotterFactory = analysisFactory->createPlotterFactory();
       plotter = plotterFactory->create();

   energ =dynamic_cast<AIDA::IDataPointSet*>(tree->find("energ")); 
   cross =dynamic_cast<AIDA::IDataPointSet*>(tree->find("cross")); 
      
   AIDA::IDataPointSet* plot_ptr;
   plot_ptr = curvFactory->createXY("curvPlot",engVal,crosVal, engErr,crosErr); 
   AIDA::IDataPointSet& curvPlot = *plot_ptr;
 

 G4cout << "............................waitting for plot ......................................." << G4endl; 
 G4cout <<" --------------------------------------------------------------------------------------" << G4endl;
 G4cout << "      curve of variation of K_CrossSections with incident energy " << G4endl; 
 G4cout << "--------------------------------------------------------------------------------------" << G4endl;

  
  plotter->show();
  plotter->currentRegion().plot(curvPlot, "cross sections");   
  
  G4cout << "Press <ENTER> to continue" << G4endl;
  G4cin.get();
  G4cin.get();

  std::cout << "Committing..." << std::endl;
  tree->commit();
  std::cout << "Closing the tree..." << std::endl;
  tree->close();

  G4cout << "END OF TEST" << G4endl;
}












