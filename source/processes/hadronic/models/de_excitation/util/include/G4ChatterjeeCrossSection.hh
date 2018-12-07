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
// V.Ivanchenko 13.04.2015
// 
// J.M. Quesada 22.04.2015 several fixes

#ifndef G4ChatterjeeCrossSection_h
#define G4ChatterjeeCrossSection_h 1

#include "globals.hh"

// index: 0-neutron, 1-proton, 2-deuteron, 3-triton, 4-He3, 5-He4
// parameters: p0, p1, p2, landa0, landa1, mu0, mu1, nu0, nu1, nu2, ra

class G4ChatterjeeCrossSection
{
public:

  static G4double ComputeCrossSection(G4double K, G4double cb, 
                                      G4double resA13, G4double amu1, 
				      G4int idx, G4int Z, G4int resA);
};

#endif
