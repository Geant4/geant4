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
// $Id: QPhysicsList.cc,v 1.2 2006-06-04 21:35:59 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QPhysicsList.cc
//
//                                         2005 Q
// ====================================================================
#include "QPhysicsList.hh"
#include "Particles.hh"
#include "PhysicsListEMstd.hh"

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////////
QPhysicsList::QPhysicsList()
  :  G4VModularPhysicsList()
////////////////////////////
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.*mm;
  SetVerboseLevel(1);

  // particles
  RegisterPhysics(new Particles);

  // EM Physics
  RegisterPhysics(new PhysicsListEMstd);
}

/////////////////////////////
QPhysicsList::~QPhysicsList()
/////////////////////////////
{
}

////////////////////////////
void QPhysicsList::SetCuts()
////////////////////////////
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}

