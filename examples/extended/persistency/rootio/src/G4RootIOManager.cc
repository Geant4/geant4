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
// $Id: G4RootIOManager.cc,v 1.4 2002-12-13 14:45:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4RootIOManager.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4RootIOManager.hh"

// Addtional Include:
#include "G4PersistencyManagerT.hh"
#include <assert.h>
#include "G4RootTransManager.hh"
#include "G4HitRootIO.hh"
#include "G4DigitRootIO.hh"
#include "G4EventRootIO.hh"
#include "G4HepMCRootIO.hh"
#include "G4MCTruthRootIO.hh"

static G4PersistencyManagerT<G4RootIOManager> rootio(G4PersistencyCenter::GetPersistencyCenter(),"ROOT");

// Implementation of Constructor #1
G4RootIOManager::G4RootIOManager(G4PersistencyCenter* pc, std::string n)
 : G4PersistencyManager(pc, n)
{
  f_G4HepMCRootIO   = new G4HepMCRootIO();
  f_G4MCTruthRootIO = G4MCTruthRootIO::GetMCTruthRootIO();
  f_G4HitRootIO     = G4HitRootIO::GetHitRootIO();
  f_G4DigitRootIO   = G4DigitRootIO::GetDigitRootIO();
  f_G4EventRootIO   = new G4EventRootIO(this);
  f_G4RootTransManager = new G4RootTransManager();
}

// Implementation of Destructor #1
G4RootIOManager::~G4RootIOManager()
{
  delete f_G4HitRootIO;
  delete f_G4DigitRootIO;
  delete f_G4EventRootIO;
  delete f_G4HepMCRootIO;
  delete f_G4MCTruthRootIO;
  delete f_G4RootTransManager;
}

// Implementation of Initialize
void G4RootIOManager::Initialize()
{
  // Initialize ROOT files for store
  f_G4RootTransManager->Initialize();

  if ( m_verbose > 2 ) {
    std::cout << " -> Root I/O is initialized." << std::endl;
  }

}

// End of G4RootIOManager.cc

