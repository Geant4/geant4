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
#ifndef G4VDNAMesh_hh
#define G4VDNAMesh_hh 1

#include "globals.hh"

class G4VDNAMesh
{
 public:
  G4VDNAMesh()          = default;
  virtual ~G4VDNAMesh() = default;
  struct Index
  {
    Index() = default;
    Index(G4int _x, G4int _y, G4int _z)
      : x(_x)
        , y(_y)
        , z(_z)
    {}
    ~Index() = default;
    G4bool operator==(const Index& rhs) const
    {
      return x == rhs.x && y == rhs.y && z == rhs.z;
    }
    G4bool operator!=(const Index& rhs) const
    {
      return x != rhs.x || y != rhs.y || z != rhs.z;
    }
    G4bool operator<(const Index& rhs) const
    {
      if(x != rhs.x)
      {
        return x < rhs.x;
      }
      else if(y != rhs.y)
      {
        return y < rhs.y;
      }
      else if(z != rhs.z)
      {
        return z < rhs.z;
      }
      else
      {
        return false;
      }
    }
    friend std::ostream& operator<<(std::ostream& s, const Index& rhs);
    G4int x = 0;
    G4int y = 0;
    G4int z = 0;
  };
  struct hashFunc
  {
    size_t operator()(const Index& k) const
    {
      size_t h1 = std::hash<G4int>()(k.x);
      size_t h2 = std::hash<G4int>()(k.y);
      size_t h3 = std::hash<G4int>()(k.z);
      return (h1 ^ (h2 << 1)) ^ h3;
    }
  };
};
#endif
