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
// class G4AuxiliaryNavServices implementation
//
// --------------------------------------------------------------------

#include "G4AuxiliaryNavServices.hh"
#include "G4GeometryTolerance.hh"

// The existence a method here allows compilers to find the inline method
// (which are the core of this class)
// [ Compilers store the inline method in the module with the first 
//   non-inline method.]
// If this method did not exist, there would be none.

void G4AuxiliaryNavServices::ReportTolerances()
{
   G4long oldPrec = G4cout.precision(16);
   
   G4cout << " Cartesian Tolerance (kCarTolerance): "
          << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()
          << " (global) " << G4endl;
   G4cout << " Radial Tolerance (kRadTolerance): "
          << G4GeometryTolerance::GetInstance()->GetRadialTolerance()
          << " (global) " << G4endl;
   G4cout << " Angular Tolerance (kAngTolerance): "
          << G4GeometryTolerance::GetInstance()->GetAngularTolerance()
          << " (global) " << G4endl;
   G4cout.precision(oldPrec); 
}
