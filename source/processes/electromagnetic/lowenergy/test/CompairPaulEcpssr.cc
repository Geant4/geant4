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
//  we try to visualize differences between Cross Sections 
// values given by exeprimental and ECPSSR methods

#include "globals.hh"
#include "G4ios.hh"
#include <vector>
#include <iostream>
#include <fstream>
#include "G4ecpssrCrossSection.hh"


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


    AIDA::IDataPointSet* ener;
    AIDA::IDataPointSet* normal;
    AIDA::IDataPointSetFactory* curvFactory;

    analysisFactory = AIDA_createAnalysisFactory();
    treeFactory = analysisFactory->createTreeFactory();
    tree = treeFactory->create("nornalisation.dat");
    curvFactory = analysisFactory->createDataPointSetFactory(*tree);

 
  
  std::vector<G4double> energy; 
  std::vector<G4double> normalisation;  
     
  std::vector<G4double> energyErr; 
  std::vector<G4double> normalisationErr; 


  G4int Z; 
  G4int zIncident; 
  G4double incidentEnergy;
  G4String fileName;
 
  G4ecpssrCrossSection* shell = new G4ecpssrCrossSection;
 
 G4cout <<"----------------------------------------------------------------------------" << G4endl; 
 G4cout <<"Enter atomic number of the incident particles (protons or alpha): " << G4endl;
 G4cout <<"----------------------------------------------------------------------------" << G4endl;
 G4cin >> zIncident;

 if (zIncident == 1)
   { fileName = "kcsPaul/kcs-";} 
  else
    {
      if (zIncident == 2)
	{ fileName = "kacsPaul/kacs-";}
	
    }

   std::vector<G4double> energies;

   for (G4double i=1.; i<3.; i=i+0.1)
     {
       energies.push_back(std::exp(i));
     } 


 G4cout <<"------------------------------------------------------------------------------" << G4endl; 
 G4cout <<"Enter atomic number of element you want to check : " << G4endl;
 G4cout <<"-------------------------------------------------------------------------------" << G4endl;
 
 G4cin >> Z;

 
 G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();

  G4VEMPDataSet* dataSet;

  dataSet = new G4EMPDataSet(Z,interpolation);
  
  dataSet->LoadData(fileName);
      

       for (size_t k=0; k<energies.size();k++)
	{
	  incidentEnergy = energies[k]/MeV;
		
          G4double logvelocity1 = log10(shell->CalculateVelocity(Z,zIncident,incidentEnergy));  	          
     	  
          G4double CS_ECPSSR = shell-> CalculateCrossSection(Z,zIncident,incidentEnergy);
         	
          G4double sigma = dataSet->FindValue(incidentEnergy) / barn;
           

	 G4cout <<"With "<<incidentEnergy/MeV <<" MeV of incident particles, K_CSection (of Z= "<<Z <<") is: "<<CS_ECPSSR<<" barn (with ECPSSR) or "<< sigma <<" barn (with PAUL_interpolt)"<< G4endl;
	 G4cout <<"----------------------------------------------------------------------------------------------------------------------------------- " << G4endl; 
        

	  G4double rapport= sigma/CS_ECPSSR;

          energy.push_back(logvelocity1);
	  normalisation.push_back(rapport);
          energyErr.push_back(0);
          normalisationErr.push_back(0);
         

         
  } 

       delete shell;
       delete dataSet;
     
      
	ener =dynamic_cast<AIDA::IDataPointSet*>(tree->find("ener")); 
	normal =dynamic_cast<AIDA::IDataPointSet*>(tree->find("normal")); 
      
        
	 
      
	  AIDA::IDataPointSet* plot_ptr;
	  plot_ptr = curvFactory->createXY("NormalisationPlot",energy,normalisation,energyErr,normalisationErr);
          AIDA::IDataPointSet& NormalisationPlot = *plot_ptr;
      

          G4cout << "........................waitting for plot of variation of (PAUL_interpolt/ECPSSR) with log(velocity)......................................." << G4endl; 
	  
          plotterFactory = analysisFactory->createPlotterFactory();
          plotter = plotterFactory->create();
          plotter->currentRegion().plot(NormalisationPlot,"Normalisation");
          plotter->show();
      
          
          G4cout << "Press <ENTER> to continue" << G4endl;
        
          G4cin.get();
	  G4cin.get();


  std::cout << "Committing..." << std::endl;
  tree->commit();
  std::cout << "Closing the tree..." << std::endl;
  tree->close();

  G4cout << "END OF TEST" << G4endl;
}

