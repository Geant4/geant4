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
// $Id: G4HepMCRootIO.cc,v 1.4 2002-12-13 14:45:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4HepMCRootIO.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4HepMCRootIO.hh"

// Implementation of Store
bool G4HepMCRootIO::Store(HepMC::GenEvent* evt)
{
  // G4RootIOManager* pm =
  //     (G4RootIOManager*) G4PersistencyCenter::GetPersistencyCenter()
  //                                           ->CurrentPersistencyManager();

  // actual implementation should come in here

  return true;
}

// Implementation of Retrieve
bool G4HepMCRootIO::Retrieve(HepMC::GenEvent*& evt, int id)
{
  // bool st = false;
  evt = 0;

  // G4RootIOManager* pm =
  //      (G4RootIOManager*) G4PersistencyCenter::GetPersistencyCenter()
  //                                            ->CurrentPersistencyManager();

  // actual implementation should come in here

  return true;
}

// End of G4HepMCRootIO.cc

