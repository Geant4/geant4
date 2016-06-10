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
// $Id: G4AllocatorList.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 

#include <iomanip>

#include "G4AllocatorList.hh"
#include "G4Allocator.hh"
#include "G4ios.hh"

G4ThreadLocal G4AllocatorList* G4AllocatorList::fAllocatorList=0;

G4AllocatorList* G4AllocatorList::GetAllocatorList()
{
  if(!fAllocatorList)
  {
    fAllocatorList = new G4AllocatorList;
  }
  return fAllocatorList;
}

G4AllocatorList* G4AllocatorList::GetAllocatorListIfExist()
{
  return fAllocatorList;
}

G4AllocatorList::G4AllocatorList()
{
}

G4AllocatorList::~G4AllocatorList()
{
  fAllocatorList = 0;
}

void G4AllocatorList::Register(G4AllocatorBase* alloc)
{
  fList.push_back(alloc);
}

void G4AllocatorList::Destroy(G4int nStat, G4int verboseLevel)
{
  std::vector<G4AllocatorBase*>::iterator itr=fList.begin();
  G4int i=0, j=0;
  G4double mem=0, tmem=0;
  if(verboseLevel>0)
  {
    G4cout << "================== Deleting memory pools ==================="
           << G4endl;
  }
  for(; itr!=fList.end();++itr)
  {
    mem = (*itr)->GetAllocatedSize();
    if(i<nStat)
    {
      i++;
      tmem += mem;
      (*itr)->ResetStorage();
      continue;
    }
    j++;
    tmem += mem;
    if(verboseLevel>1)
    {
      G4cout << "Pool ID '" << (*itr)->GetPoolType() << "', size : "
             << std::setprecision(3) << mem/1048576
             << std::setprecision(6) << " MB" << G4endl;
    }
    (*itr)->ResetStorage();
    delete *itr; 
  }
  if(verboseLevel>0)
  {
    G4cout << "Number of memory pools allocated: " << Size()
           << "; of which, static: " << i << G4endl;
    G4cout << "Dynamic pools deleted: " << j 
           << " / Total memory freed: " << std::setprecision(2)
           << tmem/1048576 << std::setprecision(6) << " MB" << G4endl;
    G4cout << "============================================================"
           << G4endl;
  }
  fList.clear();
}

G4int G4AllocatorList::Size() const
{
  return fList.size();
}
