//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4BremAngularGeneratorTest.cc,v 1.2 2003-06-16 17:00:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4BremAngularGeneratorTest.cc
//
//      Author:        Francesco Longo
// 
//      Creation date: 04 january 2001
//
//      Modifications: Luciano Pandola  (27 november 2002)
//                     Adapted in order to test G4PenelopeBremsstrahlung
//                     Minor modification in n-tuple filling
//                     Updated analysis to AIDA 3.0 
//
//                     Pedro Rodrigues/Andreia Trindade (24 March 2003)
//                     Adapted in order to test angular generators for
//                     G4LowEnergyBremsstrahlung
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"

#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4VBremAngularDistribution.hh"
#include "G4ModifiedTsai.hh"

#include "AIDA/IManagedObject.h"
//#include <memory>
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

namespace AIDA {
  class IAnalysisFactory; 
  class ITree;
  class IHistogramFactory;
  class ITupleFactory;  
  class ITuple;
  class IHistogram1D;
}

int main()
{

  // Setup

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
  AIDA::ITreeFactory* tf = af->createTreeFactory();
  AIDA::ITree* tree = tf->create("brem_angular_test.hbook","hbook",false,true);
  G4cout << "Tree store: " << tree->storeName() << G4endl;
  AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);
  AIDA::IHistogram1D* histo_1 = hf->createHistogram1D("1","Polar Angle", 360,0.,2*3.14159); 
 
  // Interactive set-up
  G4int nIterations = 1;
  G4cout << "How many interactions at generator level ?" << G4endl;
  G4cin >> nIterations;
  if (nIterations <= 0) G4Exception("Wrong input");


  G4int gType;
  G4cout << "Modified Tsai Generator [1]" << G4endl;
  G4cin >> gType;
  if ( !(gType == 1)) G4Exception("Wrong input");

  G4VBremAngularDistribution* angularDistribution = 0;

  if (gType == 1)
    {
      angularDistribution = new G4ModifiedTsai("TsaiGenerator");
    }
  angularDistribution->PrintGeneratorInformation();

  G4double initEnergy = 1*MeV; 
  G4cout << "Enter the initial electron/positron energy E (MeV)" << G4endl; 
  G4cin >> initEnergy ;
  initEnergy = initEnergy*MeV;
  
  if (initEnergy  <= 0.) G4Exception("Wrong input");

  G4double finalEnergy = 0*MeV;
  G4cout << "Enter the final electron/positron energy E (MeV)" << G4endl;
  G4cin >> finalEnergy;
  finalEnergy = finalEnergy*MeV;

  if (finalEnergy  <= 0.) G4Exception("Wrong input");

  G4int Z = 0;
  G4cout << "Enter the atomic number " << G4endl;
  G4cin >> Z;
  if(Z <= 0) G4Exception("Wrong input");

  for(G4int k = 0; k < nIterations; k++){
//    G4cout  << "Iteration number: "  <<  k << G4endl;

    G4double kineticEnergy = initEnergy;
    G4double totalEnergy = kineticEnergy + electron_mass_c2;
    G4double initial_momentum = sqrt((totalEnergy + electron_mass_c2)*kineticEnergy);

    G4double finalKineticEnergy = finalEnergy;
    G4double totalFinalEnergy = finalKineticEnergy +  electron_mass_c2;
    G4double final_momentum =  sqrt((totalFinalEnergy + electron_mass_c2)*finalKineticEnergy);

    G4double theta = angularDistribution->PolarAngle(kineticEnergy,initial_momentum,finalKineticEnergy,final_momentum,Z);
//    G4cout << theta << endl;
    histo_1->fill(theta);
  }

    G4cout << "Committing.............." << G4endl;
    tree->commit();
    G4cout << "Closing the tree........" << G4endl;
    tree->close();

    delete angularDistribution;

    delete histo_1;
    delete hf;
    delete tf;
    delete af;

   G4cout << "END OF THE MAIN PROGRAM" << G4endl;
   return 0;

}

















