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
// $Id: G4ThermalNeutrons.hh 70995 2013-06-09 00:56:34Z adotti $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4ThermalNeutrons
//
// Author: 4 May 2017 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//
// Addition of neutron thermal scattering on top of any PhysicsList

#ifndef G4ThermalNeutrons_h
#define G4ThermalNeutrons_h

#include "G4VHadronPhysics.hh"
#include "globals.hh"

class G4ThermalNeutrons : public G4VHadronPhysics {

public:
  G4ThermalNeutrons(G4int ver);
  virtual ~G4ThermalNeutrons();

  virtual void ConstructProcess();

private:
  G4int               verbose;
};

#endif






