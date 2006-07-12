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
// $Id: PhysicsListEMlowE.hh,v 1.1 2006-07-12 08:46:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   PhysicsListEMlowE.hh
//
//   Physics list for electron/positron/gamma
//   model: EM-lowE package
//
//                                         2006 Q
// ====================================================================
#ifndef PHYSICS_LIST_EMLOWE_H
#define PHYSICS_LIST_EMLOWE_H

#include "G4VUserPhysicsList.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class PhysicsListEMlowE : public G4VUserPhysicsList {

public:
  PhysicsListEMlowE();
  ~PhysicsListEMlowE();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
  virtual void SetCuts();

};

#endif

