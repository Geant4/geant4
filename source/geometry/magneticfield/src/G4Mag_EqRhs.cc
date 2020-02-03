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
// G4Mag_EqRhs implementation
//
// Created: J.Apostolakis, CERN - 13.01.1997
// --------------------------------------------------------------------

#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

const G4double G4Mag_EqRhs::fUnitConstant = 0.299792458 * (GeV/(tesla*m)); 

G4Mag_EqRhs::G4Mag_EqRhs( G4MagneticField* magField ) 
  : G4EquationOfMotion(magField)
{ 
}

G4Mag_EqRhs::~G4Mag_EqRhs()
{
}

void  
G4Mag_EqRhs::SetChargeMomentumMass( G4ChargeState particleCharge,
                                    G4double,                // MomentumXc
                                    G4double )               // particleMass
{
   G4double pcharge = particleCharge.GetCharge();
   fCof_val = pcharge*eplus*c_light ; //  B must be in Tesla
   // fCof_val = fUnitConstant*pcharge/MomentumXc; //  B must be in Tesla
   // fMass = particleMass;
}
