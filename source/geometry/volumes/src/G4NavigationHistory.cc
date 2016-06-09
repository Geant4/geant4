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
// $Id$
//
// 
// G4NavigationHistory Implementation
//
// Author: P.Kent, August 96
//
// ----------------------------------------------------------------------

#include "G4NavigationHistory.hh"
#include "G4ios.hh"

#ifndef WIN32
  // Initialise static data for the specialized memory pool
  // for the internal STL vector of histories  ...
  //
  G4ChunkIndexType* G4AllocStats::allocStat = 0;
  G4int             G4AllocStats::totSpace = 0;
  G4int             G4AllocStats::numCat = 0;
#endif

G4NavigationHistory::G4NavigationHistory()
  : fNavHistory(kHistoryMax), fStackDepth(0)
{
  Clear();
}

G4NavigationHistory::G4NavigationHistory(const G4NavigationHistory &h)
  : fNavHistory(h.fNavHistory), fStackDepth(h.fStackDepth)
{
}

G4NavigationHistory::~G4NavigationHistory()
{
}

std::ostream&
operator << (std::ostream& os, const G4NavigationHistory& nav)
{
  G4cout << "History depth=" << nav.GetDepth() << G4endl;
  for ( G4int i=0; i<=nav.GetDepth(); i++ )
  {
    os << "Level=["<<i<<"]: ";
    if( nav.GetVolume(i) != 0 )
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
