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
//
// $Id: G4PersistentRunMan.cc,v 1.10 2001/07/11 10:02:27 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// class G4PersistentRunMan 
//
// Implementation for concrete G4PersistentRunMan.
//
// History:
// 98.01.08 Y.Morita  Initial version
// 98.10.30 Y.Morita  Splitted from G4PersistencyManager

#include "G4PersistentRunMan.hh"

#include "HepODBMS/clustering/HepDbApplication.h"

#include "G4Run.hh"
#include "G4ios.hh"

G4PersistentRunMan::G4PersistentRunMan()
: f_currentRunID(0)
{;}

G4PersistentRunMan::~G4PersistentRunMan()
{;}

//----------------------------------------------------------------------------

G4bool G4PersistentRunMan::Store( HepDbApplication* dbApp,
                                  const G4Run* aRun)
{
  // Create persistent run
  f_currentPRun = new(f_container) G4PRun(aRun);

  if( f_currentPRun == 0 )
    return false;

  f_currentRunID = aRun->GetRunID();

  return true;
}

G4bool G4PersistentRunMan::Retrieve( HepDbApplication* dbApp,
                                     G4Run*& aRun)
{
  G4bool theStatus = false;
  aRun = 0;

  ooItr(G4PRun) pRun_iterator;

  // set the new scan scope if the DB/Container has been changed by
  // G4PersistencyManager
  if( f_container != f_currentContainer )
  {
    f_currentContainer = f_container;
    pRun_iterator.scan(f_currentContainer);
  }

  // Retrieve "next" G4PRun in this scope and make a G4Run object
  if( pRun_iterator.next() )
  {
    G4Run* run = pRun_iterator->MakeTransientObject();
    if( run != 0 )
    {
      aRun = run;
      theStatus = true;
      f_currentRunID = aRun->GetRunID();
    }
  }

  return theStatus;
}

