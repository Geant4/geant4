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
// $Id: G4PhotoElectricAngularGeneratorTest.cc,v 1.1 2004-05-12 09:24:53 silvarod Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4PhotoElectricAngularGeneratorTest.cc
//
//      Author:        Francesco Longo
// 
//      Creation date: 04 january 2001
//
//      Modifications: Luciano Pandola  (27 november 2002)
//                     Adapted in order to test G4PenelopePhotoElectricsstrahlung
//                     Minor modification in n-tuple filling
//                     Updated analysis to AIDA 3.0 
//                     Pedro Rodrigues (10 May 2004)
//                     Adapted in order to test angular generators for
//                     G4LowEnergyPhotoElectricsstrahlung
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"

#include "fstream"
#include "iomanip"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4VPhotoElectricAngularDistribution.hh"
#include "G4PhotoElectricAngularGenerator462.hh"
#include "G4PhotoElectricAngularGeneratorStandard.hh"


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

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  G4std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );
  G4std::auto_ptr< AIDA::ITreeFactory > tf (af->createTreeFactory());
  G4std::auto_ptr< AIDA::ITree > tree (tf->create("brem_angular_test.hbook","hbook",false,true));
  G4cout << "Tree store: " << tree->storeName() << G4endl;
  G4std::auto_ptr< AIDA::IHistogramFactory > hf (af->createHistogramFactory(*tree));
  G4std::auto_ptr< AIDA::IHistogram1D> histo_1 (hf->createHistogram1D("1","Cos Polar Angle", 100,0.,1.0)); 
  G4std::auto_ptr< AIDA::IHistogram1D> histo_2 (hf->createHistogram1D("2","Azimuthal Angle", 100,0.,3.14159)); 

  // Interactive set-up
  G4cout << "How many interactions at generator level ?" << G4endl;
  G4cin >> nIterations;
  if (nIterations <= 0) G4Exception("Wrong input");


  G4int gType = 0;
  G4cout << "Toy Model [1]" << G4endl;
  G4cout << "Standard  [2]" << G4endl;
  G4cin >> gType;
  if ( !(gType < 3)) G4Exception("Wrong input");

  G4VPhotoElectricAngularDistribution* angularDistribution = 0;

  if (gType == 1)
    {
      angularDistribution = new G4PhotoElectricAngularGenerator462("PhotoElectricAngularGenerator462");
    }

  if(gType == 2)
    {
      angularDistribution = new G4PhotoElectricAngularGeneratorStandard("PhotoElectricAngularGeneratorStandard");
    }

  angularDistribution->PrintGeneratorInformation();

  G4ThreeVector direction(0,0,1); //photon incident energy
  G4double initEnergy = 0;
  G4cout << "Enter the photon energy E (keV)" << G4endl; 
  G4cin >> initEnergy ;
  initEnergy = initEnergy*keV;
  
  if (initEnergy  <= 0.) G4Exception("Wrong input");

  for(G4int i = 0; i < nIterations; i++){
    G4ThreeVector photondirection = angularDistribution->GetPhotoElectronDirection(direction, initEnergy);
    histo->Fill((photondirection.getZ()));
    G4double h = sqrt(1-photondirection.getZ()*photondirection.getZ());
    if(h != 0) histo2->Fill(acos(photondirection.getX()/h));

    G4cout << "Committing.............." << G4endl;
    tree->commit();
    G4cout << "Closing the tree........" << G4endl;
    tree->close();
  }
    delete angularDistribution;

   G4cout << "END OF THE MAIN PROGRAM" << G4endl;
   return 0;

}

















