//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: PhysicsListEMstd.hh,v 1.3 2006-06-04 21:36:35 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   PhysicsListEMstd.hh
//
//   Physics list for electron/positron/gamma
//   EM-standard package
//
//                                         2006 Q
// ====================================================================
#ifndef PHYSICS_LIST_EMSTD_H
#define PHYSICS_LIST_EMSTD_H

#include "G4VUserPhysicsList.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class PhysicsListEMstd: public G4VUserPhysicsList {
public:
  PhysicsListEMstd();
  ~PhysicsListEMstd();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
  virtual void SetCuts();

};

#endif
