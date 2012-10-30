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
//      File name:     G4PiMinusStopCo
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
// -------------------------------------------------------------------

#include "G4ios.hh"

#include "G4PiMinusStopCo.hh"

#include <vector>

#include "globals.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ParticleTypes.hh"
#include "G4ReactionKinematics.hh"
#include "G4DynamicParticleVector.hh"
#include "G4LorentzVector.hh"
#include "G4PiMinusStopMaterial.hh"
#include "G4DistributionGenerator.hh"

// np/pp production ratio
// Experimental values: 
// R(np/pp) (E. Gadioli et al., Phys Rev C36 (1987) 741
G4double G4PiMinusStopCo::npRatio = 2.5;


 
// Average numbers of final nucleons detected, for N-pair absorption
// (Hartmann et al., Nucl. Phys. A300 (1978) 345
G4double G4PiMinusStopCo::nFinalNucleons = 1.38;

// Kinetic energy (MeV) distributions measured for coincident nucleon 
// emission
// P. Heusi et al., Nucl. Phys. A407(1983) 429

G4int G4PiMinusStopCo::eKinEntries = 11;

G4double G4PiMinusStopCo::eKinData[11] = {0.085, 0.09, 0.09, 0.15, 
				      0.1, 0.09, 0.08, 
				      0.04,  0.03,  0.02, 0.01};

G4double G4PiMinusStopCo::eKin[12] = { 15., 17.5, 25.,  33.,  
				   42.,  52.,  62., 
				   75., 85., 95., 105. };


// Opening angle distributions measured for coincident nucleon emission
// (P.Heusi et al., Nucl. Phys. A407 (1983) 429

G4int G4PiMinusStopCo::angleEntries = 7;

G4double G4PiMinusStopCo::angleData[7] = 
{6., 8., 9., 10., 25., 40., 45. };

G4double G4PiMinusStopCo::angle[8] = { 0.5, 0.7, 0.87, 1.4, 2.1, 
				   2.44, 2.8, 3.1415927 };



// Constructor

G4PiMinusStopCo::G4PiMinusStopCo()
  
{
  // Cluster size: nucleon pair, alpha, triton etc.
  // First implementation: interaction with nucleon pair only
  _clusterSize = 2;

  // R ratio
  theR = 1. / (1. + npRatio);

  _definitions = new std::vector<G4ParticleDefinition*>();
  _momenta = new std::vector<G4LorentzVector*>();

  std::vector<double> eKinVector;
  std::vector<double> eKinDataVector;
  int i;
  for (i=0; i<eKinEntries; i++)
    {
      eKinVector.push_back(eKin[i]);
      eKinDataVector.push_back(eKinData[i]);
    }
  eKinVector.push_back(eKin[eKinEntries]);
  _distributionE = new G4DistributionGenerator(eKinVector,eKinDataVector);

  std::vector<double> angleVector;
  std::vector<double> angleDataVector;
  for (i=0; i<angleEntries; i++)
    {
      angleVector.push_back(angle[i]);
      angleDataVector.push_back(angleData[i]);
    }
  angleVector.push_back(angle[angleEntries]);
  _distributionAngle = new G4DistributionGenerator(angleVector,angleDataVector);
}


// Destructor

G4PiMinusStopCo::~G4PiMinusStopCo()
{}

G4double G4PiMinusStopCo::FinalNucleons()
{
  return nFinalNucleons;
}

