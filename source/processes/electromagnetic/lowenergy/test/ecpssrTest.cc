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

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4ecpssrCrossSection.hh"

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
 

  G4int zIncident;
  G4int zTarget; 
  G4double energyIncident;
 
  G4ecpssrCrossSection* sigma_ecpssr = new G4ecpssrCrossSection;

 

  std::vector<G4double> energies; 
  for (G4double i=-3.; i<3.0; i=i+0.25)
     {
       energies.push_back(std::exp(i));
     } 
   
 G4cout <<""<<G4endl;
 G4cout <<"----------------------------------------------------------------------------" << G4endl; 
 G4cout <<"Enter atomic number of the incident particles (protons or alpha): " << G4endl;
 G4cout <<"----------------------------------------------------------------------------" << G4endl;
 G4cin >> zIncident;


 G4cout <<""<<G4endl;
 G4cout <<"----------------------------------------------------------------------------" << G4endl; 
 G4cout <<"Enter atomic number of element you want to check its K_Shell CrossSection: " << G4endl;
 G4cout <<"----------------------------------------------------------------------------" << G4endl;
 G4cin >> zTarget;



 
      for (size_t k=0; k<energies.size();k++)                       //Cross section for each incident energy
	{

	  energyIncident = energies[k]/MeV;
           
          
         
	  G4double K_CrossSection = sigma_ecpssr->CalculateCrossSection(zTarget,zIncident,energyIncident);
        
	    
	 G4cout <<"at incident Particles of "<<energyIncident/MeV <<" MeV, the K_CrossSection, of the element Z=" <<zTarget <<", is "<<K_CrossSection<<" barn" << G4endl;
	 G4cout <<"-----------------------------------------------------------------------------------------" << G4endl; 
         
          engVal.push_back(energyIncident);
          crosVal.push_back(K_CrossSection);
          engErr.push_back(0);
          crosErr.push_back(0);
	  
	   }   
	  
   

       delete sigma_ecpssr;
       
        
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
