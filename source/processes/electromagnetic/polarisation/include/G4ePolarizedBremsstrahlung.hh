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
// $Id: G4ePolarizedBremsstrahlung.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ePolarizedBremsstrahlung
//
// Author:        Karim Laihem based on code by Michel Maire
//
// Creation date: 01.05.2005
//
// Modifications:
// 21-08-06 Modified to work in g4.8.1 framework (A.Schaelicke)
//
// Class Description:
//
// polarized version of G4eBremsstrahlung

#ifndef G4ePolarizedBremsstrahlung_h
#define G4ePolarizedBremsstrahlung_h 1

#include "G4eBremsstrahlung.hh"

class G4ePolarizedBremsstrahlung : public G4eBremsstrahlung
{

public:

  explicit G4ePolarizedBremsstrahlung(const G4String& name = "pol-eBrem");

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*) override;

};
#endif
