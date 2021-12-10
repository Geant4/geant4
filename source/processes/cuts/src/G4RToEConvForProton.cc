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
// G4RToEConvForProton class implementation
//
// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------

#include "G4RToEConvForProton.hh"
#include "G4ParticleTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// --------------------------------------------------------------------
G4RToEConvForProton::G4RToEConvForProton() 
  : G4VRangeToEnergyConverter()
{    
  theParticle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  if (theParticle == nullptr)
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4RToEConvForProton::G4RToEConvForProton() - ";
      G4cout << "Proton is not defined !!" << G4endl;
    }
#endif
  }
  else 
  {
    fPDG = theParticle->GetPDGEncoding();
  }
}

// --------------------------------------------------------------------
G4RToEConvForProton::~G4RToEConvForProton()
{}

// --------------------------------------------------------------------
G4double G4RToEConvForProton::Convert(const G4double rangeCut, 
                                      const G4Material* )
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>3)
  {
    G4cout << "G4RToEConvForProton::Convert() - ";
    G4cout << " with Range Cut " << rangeCut/mm << "[mm]" << G4endl;
  }
#endif
  // Simple formula - range = Ekin/(100*keV)*(1*mm);
  return (rangeCut/(1.0*CLHEP::mm)) * (100.0*CLHEP::keV); 
}

// --------------------------------------------------------------------
G4double G4RToEConvForProton::ComputeValue(const G4int, const G4double)
{
  return 0.0;
}

// --------------------------------------------------------------------
