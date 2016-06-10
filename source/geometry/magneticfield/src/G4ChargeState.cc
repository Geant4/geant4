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
// History
// - First version: Apr 10, 2013  John Apostolakis, Peter Gumplinger
// - Modified:
//
//
// $Id: G4FieldTrack.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// -------------------------------------------------------------------

#include "G4ChargeState.hh"

void G4ChargeState::SetChargeSpinMoments(G4double charge,
                                         G4double spin,
                                         G4double magnetic_dipole_moment,
                                         G4double electric_dipole_moment,
                                         G4double magnetic_charge )
   //  Revise the charge and potentially all moments.
   //   By default do not change mdm, edm, mag charge.
{
   fCharge = charge;
   fSpin   = spin;
   if( magnetic_dipole_moment < DBL_MAX) fMagn_dipole= magnetic_dipole_moment;
   if( electric_dipole_moment < DBL_MAX) fElec_dipole= electric_dipole_moment;
   if( magnetic_charge < DBL_MAX)        fMagneticCharge= magnetic_charge;
}
