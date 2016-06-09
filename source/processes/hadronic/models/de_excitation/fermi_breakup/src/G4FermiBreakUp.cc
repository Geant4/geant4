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
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4FermiBreakUp.hh"
#include "G4HadronicException.hh"

G4FermiBreakUp::G4FermiBreakUp()
{
}

G4FermiBreakUp::G4FermiBreakUp(const G4FermiBreakUp &) : G4VFermiBreakUp()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4FermiBreakUp::copy_constructor meant to not be accessable");
}


G4FermiBreakUp::~G4FermiBreakUp()
{
}


const G4FermiBreakUp & G4FermiBreakUp::operator=(const G4FermiBreakUp &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4FermiBreakUp::operator= meant to not be accessable");
  return *this;
}


G4bool G4FermiBreakUp::operator==(const G4FermiBreakUp &) const
{
  return false;
}

G4bool G4FermiBreakUp::operator!=(const G4FermiBreakUp &) const
{
  return true;
}



G4FragmentVector * G4FermiBreakUp::BreakItUp(const G4Fragment &theNucleus)
{
  // CHECK that Excitation Energy > 0
  if (theNucleus.GetExcitationEnergy() <= 0) 
    {
      G4FragmentVector * theResult = new G4FragmentVector;
      theResult->push_back(new G4Fragment(theNucleus));
      return theResult;
    }
  
  // Total energy of nucleus in nucleus rest frame 
  G4double TotalEnergyRF = theNucleus.GetMomentum().m();
  //   G4double TotalEnergyRF = theNucleus.GetExcitationEnergy() +
  //     G4ParticleTable::GetParticleTable()->GetIonTable()->
  //     GetIonMass(static_cast<G4int>(theNucleus.GetZ()),static_cast<G4int>(theNucleus.GetA()));
  
  
  G4FermiConfigurationList theConfigurationList;
  
  
  // Split the nucleus
  G4bool Split = theConfigurationList.Initialize(static_cast<G4int>(theNucleus.GetA()), 
						 static_cast<G4int>(theNucleus.GetZ()),
						 TotalEnergyRF);
  if ( !Split ) 
    {
      G4FragmentVector * theResult = new G4FragmentVector;
      theResult->push_back(new G4Fragment(theNucleus));
      
      return theResult;
    }

  // Chose a configuration
  G4FermiConfiguration theConfiguration(theConfigurationList.ChooseConfiguration());
  
  
  // Get the fragments corresponding to chosen configuration.
  G4FragmentVector * theResult = theConfiguration.GetFragments(theNucleus);
#ifdef PRECOMPOUND_TEST
  for (G4FragmentVector::iterator i = theResult->begin(); i != theResult->end(); i++)
    {
      (*i)->SetCreatorModel("G4FermiBreakUp");
    }
#endif
  return theResult;
  
}


