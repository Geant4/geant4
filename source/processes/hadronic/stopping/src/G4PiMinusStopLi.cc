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
//      File name:     G4PiMinusStopLi
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
// -------------------------------------------------------------------

#include "G4ios.hh"

#include "G4PiMinusStopLi.hh"

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
// R(np/pp) = 6.3 +- 1.4 (P.Heusi et al., Nucl. Phys. A407 (1983) 429
G4double G4PiMinusStopLi::npRatio = 6.3;
 
// Average numbers of final nucleons detected, for N-pair absorption
// (P.Heusi et al., Nucl. Phys. A407 (1983) 429
G4double G4PiMinusStopLi::nFinalNucleons = 1.9;

// Kinetic energy (MeV) distributions measured for coincident nucleon 
// emission
// P. Heusi et al., Nucl. Phys. A407(1983) 429

G4int G4PiMinusStopLi::eKinEntries = 21;

G4double G4PiMinusStopLi::eKinData[21] = { 0.0018, 0.0025, 0.003, 0.045, 
				      0.007, 0.014, 0.023, 
				      0.4,  0.09,  0.18, 
				      0.25,  0.3,  0.25, 0.2, 
				      0.18,  0.08,  0.05, 
				      0.023, 0.012,  0.007, 0.02};

G4double G4PiMinusStopLi::eKin[22] = { 15., 17.5, 22.5,  27.5,  
				   32.5,  37.5,  42.5,  
				   47.5, 52.5,  57.5,  
				   62.5, 67.5,  72.5,  77.5, 
				   82.5, 87.5,  92.5,  
				   97.5, 102.5, 105. };


// Opening angle distributions measured for coincident nucleon emission
// (P.Heusi et al., Nucl. Phys. A407 (1983) 429

G4int G4PiMinusStopLi::angleEntries = 7;

G4double G4PiMinusStopLi::angleData[7] = 
{ 0.17, 0.4, 0.7, 1.1, 1.3, 20., 70. };

G4double G4PiMinusStopLi::angle[8] = { 1.308997, 1.570796, 1.832596, 2.094395, 
				  2.356194, 2.617994, 2.967060, 3.1415927 };



// Constructor

G4PiMinusStopLi::G4PiMinusStopLi()
  
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

G4PiMinusStopLi::~G4PiMinusStopLi()
{}

G4double G4PiMinusStopLi::FinalNucleons()
{
  return nFinalNucleons;
}

