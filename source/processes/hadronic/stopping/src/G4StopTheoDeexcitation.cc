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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StopTheoDeexcitation.cc,v 1.8 2002-06-07 15:30:58 jwellisc Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
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
#include "G4NucleiPropertiesTable.hh"
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
  // theHandler.SetMinEForMultiFrag(3*MeV);
  theHandler.SetMinEForMultiFrag(300*GeV);

  // Deexcite the nucleus 

  G4double atomicMass = G4NucleiPropertiesTable::GetAtomicMass(static_cast<G4int>(Z),static_cast<G4int>(A));
  G4double m = atomicMass + excitation;
  G4double pMag = p.mag();
  G4LorentzVector initialMomentum(p.x(),p.y(),p.z(),sqrt(pMag*pMag + m*m));
  G4Fragment theExcitedNucleus(static_cast<G4int>(A),static_cast<G4int>(Z),initialMomentum);

  //  theExcitedNucleus.SetA(A);
  //  theExcitedNucleus.SetZ(Z);
  //  theExcitedNucleus.SetExcitationEnergy(excitation); 
  //  theExcitedNucleus.SetMomentum(initialMomentum);

  //  G4cout << "Theo input " << A << " " << Z << " " 
  //	 << pMag << " " << atomicMass << G4endl
  //	 << "Theo -     " << excitation << " " << initialMomentum.mag() << G4endl
  //	 << "Fragment - " << theExcitedNucleus.GetExcitationEnergy() << " "
  //	 << theExcitedNucleus.GetMomentum().mag() << G4endl;

  return theHandler.BreakItUp(theExcitedNucleus);
}


