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
// G4NavigationHistory Implementation
//
// Author: P.Kent, August 96
//
// ----------------------------------------------------------------------

#include "G4ios.hh"
#include "G4NavigationHistory.hh"

G4Allocator<G4NavigationHistory>*& aNavigHistoryAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4NavigationHistory>* _instance = nullptr;
  return _instance;
}

G4NavigationHistory::G4NavigationHistory()
  : fStackDepth(0)
{
  fNavHistory = G4NavigationHistoryPool::GetInstance()->GetLevels();
  Clear();
}

G4NavigationHistory::G4NavigationHistory(const G4NavigationHistory &h)
{
  fNavHistory = G4NavigationHistoryPool::GetInstance()->GetLevels();
  if( GetMaxDepth() != h.GetMaxDepth() )
  {
    fNavHistory->resize( h.GetMaxDepth() );
  }

  for ( G4long ilev=G4long(h.fStackDepth); ilev>=0; --ilev )
  { 
    (*fNavHistory)[ilev] = (*h.fNavHistory)[ilev];
  }
  fStackDepth = h.fStackDepth;
}

G4NavigationHistory::~G4NavigationHistory()
{
  G4NavigationHistoryPool::GetInstance()->DeRegister(fNavHistory);
}

std::ostream&
operator << (std::ostream& os, const G4NavigationHistory& nav)
{
  os << "History depth=" << nav.GetDepth() << G4endl;
  for ( G4int i=0; i<=(G4int)nav.GetDepth(); ++i )
  {
    os << "Level=["<<i<<"]: ";
    if( nav.GetVolume(i) != nullptr )
    {
      os << "Phys Name=["<< nav.GetVolume(i)->GetName()
         << "] Type=[";
      switch(nav.GetVolumeType(i))
      {
        case kNormal:
          os << "N";
          break;
        case kReplica:
          os << "R" << nav.GetReplicaNo(i);
          break;
        case kParameterised:
          os << "P" << nav.GetReplicaNo(i);
          break;
        case kExternal:
          os << "E" << nav.GetReplicaNo(i);
          break;
      }
      os << "]";
    }
    else
    {
      os << "Phys = <Null>";
    }
    os << G4endl;
  }
  return os;
}
