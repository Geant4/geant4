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
// $Id: G4BremAngularGeneratorTest.cc,v 1.5 2004-03-14 21:57:25 silvarod Exp $
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
//                     Pedro Rodrigues (24 March 2003)
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
#include "G4Generator2BS.hh"
#include "G4Generator2BN.hh"

#include "AIDA/IManagedObject.h"
#include <memory>
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

int main()
{

  // Setup

  // -------------------------------------------------------------------

  // ---- HBOOK initialization
  std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );
  std::auto_ptr< AIDA::ITreeFactory > tf (af->createTreeFactory());
  std::auto_ptr< AIDA::ITree > tree (tf->create("brem_angular_test.hbook","hbook",false,true));
  G4cout << "Tree store: " << tree->storeName() << G4endl;
  std::auto_ptr< AIDA::IHistogramFactory > hf (af->createHistogramFactory(*tree));
  std::auto_ptr< AIDA::IHistogram1D> histo_1 (hf->createHistogram1D("1","Polar Angle", 100,0.,3.14159)); 
  // Interactive set-up
  G4int nIterations = 1;
  G4cout << "How many interactions at generator level ?" << G4endl;
  G4cin >> nIterations;
  if (nIterations <= 0) G4Exception("Wrong input");


  G4int gType;
  G4cout << "Modified Tsai Generator [1]" << G4endl;
  G4cout << "2BS Generator [2]" << G4endl;
  G4cout << "2BN Generator [3]" << G4endl;
  G4cin >> gType;
  if ( !(gType < 4)) G4Exception("Wrong input");

  G4VBremAngularDistribution* angularDistribution = 0;

  if (gType == 1)
    {
      angularDistribution = new G4ModifiedTsai("TsaiGenerator");
    }

  if(gType == 2)
    {
      angularDistribution = new G4Generator2BS("2BSGenerator");
    }

  if(gType == 3)
    {
      angularDistribution = new G4Generator2BN("2BNGenerator");
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

  if (finalEnergy  < 0.) G4Exception("Wrong input");

  G4int Z = 0;
  G4cout << "Enter the atomic number " << G4endl;
  G4cin >> Z;
  if(Z <= 0) G4Exception("Wrong input");

  G4double kineticEnergy = initEnergy;
  G4double finalKineticEnergy = finalEnergy;

  for(G4int k = 0; k < nIterations; k++){
    G4double theta = angularDistribution->PolarAngle(kineticEnergy,finalKineticEnergy,Z);
    histo_1->fill(theta);
  }

    G4cout << "Committing.............." << G4endl;
    tree->commit();
    G4cout << "Closing the tree........" << G4endl;
    tree->close();

    delete angularDistribution;

   G4cout << "END OF THE MAIN PROGRAM" << G4endl;
   return 0;

}

















