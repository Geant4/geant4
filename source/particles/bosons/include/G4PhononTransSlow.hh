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
/// \file particles/phonons/include/G4PhononTransSlow.hh
/// \brief Definition of the G4PhononTransSlow class
//
// $Id: G4PhononTransSlow.hh 75122 2013-10-28 09:51:40Z gcosmo $
//

#ifndef G4PhononTransSlow_h
#define G4PhononTransSlow_h 1

#include "G4ParticleDefinition.hh"

class G4PhononTransSlow : public G4ParticleDefinition {
private:
  static G4PhononTransSlow* theInstance;

private:
  G4PhononTransSlow() {;}

public:
  virtual ~G4PhononTransSlow () {;}
  
  static G4PhononTransSlow* Definition();
  static G4PhononTransSlow* PhononDefinition();
};

#endif

