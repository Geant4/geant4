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
//      File name:     G4StopTheoDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
//      Modifications: 
//      MGP            4 July 1998 Modified parameters to force evaporation
//                     
// -------------------------------------------------------------------

#include "G4ios.hh"

#include "G4StopTheoDeexcitation.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4NucleiProperties.hh"
#include "G4Fragment.hh"
#include "G4ExcitationHandler.hh"
#include "G4DynamicParticleVector.hh"

// Constructor

G4StopTheoDeexcitation::G4StopTheoDeexcitation()
{}

// Destructor

G4StopTheoDeexcitation::~G4StopTheoDeexcitation()
{}

G4ReactionProductVector* G4StopTheoDeexcitation::BreakUp(G4double A, G4double Z, 
							 G4double excitation, 
							 const G4ThreeVector& p)
{

  G4ExcitationHandler theHandler;

  // MF and FB parameters modified by MGP to force evaporation 
  // Max A and Z values for use Fermi Breakup
  // theHandler.SetMaxAandZForFermiBreakUp(16, 10);
  //  theHandler.SetMaxAandZForFermiBreakUp(2, 1);
  // Min excitation energy (per nucleon) for use MultiFrag

  theHandler.SetMinEForMultiFrag(300*GeV);

  // Deexcite the nucleus 

  G4double atomicMass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(A),static_cast<G4int>(Z));
  G4double m = atomicMass + excitation;
  G4double pMag = p.mag();
  G4LorentzVector initialMomentum(p.x(),p.y(),p.z(),std::sqrt(pMag*pMag + m*m));
  G4Fragment theExcitedNucleus(static_cast<G4int>(A),static_cast<G4int>(Z),initialMomentum);

  return theHandler.BreakItUp(theExcitedNucleus);
}


