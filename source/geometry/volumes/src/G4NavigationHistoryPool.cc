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
// $Id:$
//
// G4NavigationHistoryPool
//
// Implementation for singleton container
//
// History:
// 07.05.14 G.Cosmo Initial version
// --------------------------------------------------------------------

#include "globals.hh"
#include "G4NavigationHistoryPool.hh"

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4ThreadLocal G4NavigationHistoryPool* G4NavigationHistoryPool::fgInstance = 0;

// ***************************************************************************
// Private constructor: Construct underlying containers
// ***************************************************************************
//
G4NavigationHistoryPool::G4NavigationHistoryPool()
{
  fPool.reserve(512);
  fFree.reserve(512);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4NavigationHistoryPool::~G4NavigationHistoryPool() 
{
  Clean(); fgInstance = 0;
}

// ***************************************************************************
// Delete all elements from the pool
// ***************************************************************************
//
void G4NavigationHistoryPool::Clean()
{
  for(size_t i=0; i<fPool.size(); ++i)
  {
    delete fPool[i];
  }
  fPool.clear();
  fFree.clear();
}

// ***************************************************************************
// Print number of entries
// ***************************************************************************
//
void G4NavigationHistoryPool::Print() const
{
#ifdef G4VERBOSE
  G4cout << "Total navigation history collections cleaned: "
         << fPool.size() << G4endl;
#endif
}

// ***************************************************************************
// Delete all elements from the pool
// ***************************************************************************
//
void G4NavigationHistoryPool::Reset()
{
  for(size_t i=0; i<fPool.size(); ++i)
  {
    fPool[i] = 0;
  }
  for(size_t j=0; j<fFree.size(); ++j)
  {
    fFree[j] = 0;
  }
}

// ***************************************************************************
// Return ptr to Store, setting if necessary
// ***************************************************************************
//
G4NavigationHistoryPool* G4NavigationHistoryPool::GetInstance()
{
  if (!fgInstance)
  {
    fgInstance = new G4NavigationHistoryPool;
  }
  return fgInstance;
}
