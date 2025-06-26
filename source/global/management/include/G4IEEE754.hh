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
#ifndef G4IEEE754_hh
#define G4IEEE754_hh 1

#include "G4Types.hh"

#include <cstdint>

namespace G4IEEE754
{
//----------------------------------------------------------------------------
// Used to switch between different type of interpretations of the data
// (64 bits)
//
union ieee754
{
    ieee754() = default;
    ieee754(G4double thed) { d = thed; };
    ieee754(uint64_t thell) { ll = thell; };
    ieee754(G4float thef) { f[0] = thef; };
    ieee754(uint32_t thei) { i[0] = thei; };
    G4double d;
    G4float f[2];
    uint32_t i[2];
    uint64_t ll;
    uint16_t s[4];
};

//----------------------------------------------------------------------------
// Converts a double to an unsigned long long
//
inline uint64_t dp2uint64(G4double x)
{
  ieee754 tmp;
  tmp.d = x;
  return tmp.ll;
}

//----------------------------------------------------------------------------
// Converts an unsigned long long to a double
//
inline G4double uint642dp(uint64_t ll)
{
  ieee754 tmp;
  tmp.ll = ll;
  return tmp.d;
}

//----------------------------------------------------------------------------
// Converts an int to a float
//
inline G4float uint322sp(G4int x)
{
  ieee754 tmp;
  tmp.i[0] = x;
  return tmp.f[0];
}

//----------------------------------------------------------------------------
// Converts a float to an int
//
inline uint32_t sp2uint32(G4float x)
{
  ieee754 tmp;
  tmp.f[0] = x;
  return tmp.i[0];
}
}  // namespace G4IEEE754

#endif
