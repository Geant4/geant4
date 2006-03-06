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
// $Id: G4PhotoElectricAngularGeneratorTest.cc,v 1.3 2006-03-06 16:41:11 silvarod Exp $
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

#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

//#include "fstream"
//#include "iomanip"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4VPhotoElectricAngularDistribution.hh"
#include "G4PhotoElectricAngularGenerator462.hh"
#include "G4PhotoElectricAngularGeneratorStandard.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"

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

G4ThreeVector SetPerpendicularVector(G4ThreeVector& a)
{
  G4double dx = a.x();
  G4double dy = a.y();
  G4double dz = a.z();
  G4double x = dx < 0.0 ? -dx : dx;
  G4double y = dy < 0.0 ? -dy : dy;
  G4double z = dz < 0.0 ? -dz : dz;
  if (x < y) {
    return x < z ? G4ThreeVector(-dy,dx,0) : G4ThreeVector(0,-dz,dy);
  }else{
    return y < z ? G4ThreeVector(dz,0,-dx) : G4ThreeVector(-dy,dx,0);
  }
}

int main()
{

  // -------------------------------------------------------------------

  // ---- HBOOK initialization

  AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
  AIDA::ITreeFactory* tf = af->createTreeFactory();
  AIDA::ITree* tree = tf->create("phot_angular_test.hbook","hbook",false,true);
  G4cout << "Tree store: " << tree->storeName() << G4endl;
  AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);
  AIDA::IHistogram1D* histo_1 = hf->createHistogram1D("1","Cos Polar Angle", 100,0.,3.14159); 
  AIDA::IHistogram1D* histo_2 = hf->createHistogram1D("2","Azimuthal Angle", 100,0.,3.14159); 
  AIDA::IHistogram2D* histo_3 = hf->createHistogram2D("3","2d", 100,0.,3.14159,200,-3.14159,3.14159); 

  // Interactive set-up
  G4cout << "How many interactions at generator level ?" << G4endl;
  G4int nIterations = 0;
  G4cin >> nIterations;
  if (nIterations <= 0) G4Exception("Wrong input");

  G4int gType = 0;
  G4cout << "Toy Model [1]" << G4endl;
  G4cout << "Standard  [2]" << G4endl;
  G4cout << "Polarized [3]" << G4endl;
  G4cin >> gType;
  if ( !(gType < 4)) G4Exception("Wrong input");

  G4VPhotoElectricAngularDistribution* angularDistribution = 0;

  if (gType == 1)
    {
      angularDistribution = new G4PhotoElectricAngularGenerator462("PhotoElectricAngularGenerator462");
    }

  if(gType == 2)
    {
      angularDistribution = new G4PhotoElectricAngularGeneratorStandard("PhotoElectricAngularGeneratorStandard");
    }

  if(gType == 3)
    {
      angularDistribution = new G4PhotoElectricAngularGeneratorPolarized("PhotoElectricAngularGeneratorPolarized");
    }

  angularDistribution->PrintGeneratorInformation();

  G4ThreeVector direction(0,0,1); //photon incident energy
  G4double initEnergy = 0;
  G4cout << "Enter the photon energy E (keV)" << G4endl; 
  G4cin >> initEnergy ;
  initEnergy = initEnergy*keV;
  
  if (initEnergy  <= 0.) G4Exception("Wrong input");
  G4cout << "Shell level" << G4endl;
  G4int level;
  G4cin >> level;

  for(G4int i = 0; i < nIterations; i++){
    
    // random poll
    G4ThreeVector d0 = direction.unit();
    G4ThreeVector a1 = SetPerpendicularVector(d0); //different orthogonal
    G4ThreeVector a0 = a1.unit(); // unit vector
//    G4double rand1 = G4UniformRand();
    G4double rand1 = 0;
    G4double angle = twopi*rand1; // random polar angle
    G4ThreeVector b0 = d0.cross(a0); // cross product
    G4ThreeVector c;
    c.setX(std::cos(angle)*(a0.x())+std::sin(angle)*b0.x());
    c.setY(std::cos(angle)*(a0.y())+std::sin(angle)*b0.y());
    c.setZ(std::cos(angle)*(a0.z())+std::sin(angle)*b0.z());
    G4ThreeVector pol = c.unit();

    G4ThreeVector photondirection = angularDistribution->GetPhotoElectronDirection(direction, initEnergy, pol, level);


    histo_1->fill(photondirection.theta());
    histo_2->fill(photondirection.phi());
    histo_3->fill(photondirection.theta(),photondirection.phi());

  }

    G4cout << "Committing.............." << G4endl;
    tree->commit();
    G4cout << "Closing the tree........" << G4endl;
    tree->close();
  
   delete angularDistribution;
   G4cout << "END OF THE MAIN PROGRAM" << G4endl;
   return 0;

}

















