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
#ifndef MOLECULAR_UTILITY_FUNCTIONS_HH
#define MOLECULAR_UTILITY_FUNCTIONS_HH

#include "globals.hh"
#include <string>
#include <vector>
#include <array>
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace utility
{
  std::vector<G4String>& Split(const G4String&, char, std::vector<G4String>&);

  std::vector<G4String> Split(const G4String&, char);

  G4String Get_seperated_element(const G4String&, char, G4int);

  std::array<G4String, 4> Get_four_elements(const G4String&, char);

  G4bool Path_exists(const G4String&);

  G4String Getcwd();

  G4double Min(const G4ThreeVector&);

  // sign function
  template <typename T>
  [[maybe_unused]] inline constexpr G4int Signum(T, std::false_type);

  template <typename T>
  inline constexpr int Signum(T, std::true_type);

  template <typename T>
  inline constexpr int Signum(T);

  // find last (highest) significant bit
  // linux implementation, copied for portability
  // http://lxr.linux.no/#linux+v2.6.26.5/include/asm-generic/bitops/fls.h
  static inline int Fls(int x)
  {
    int r = 32;

    if(!x)
      return 0;
    if(!(x & 0xffff0000u))
    {
      x <<= 16;
      r -= 16;
    }
    if(!(x & 0xff000000u))
    {
      x <<= 8;
      r -= 8;
    }
    if(!(x & 0xf0000000u))
    {
      x <<= 4;
      r -= 4;
    }
    if(!(x & 0xc0000000u))
    {
      x <<= 2;
      r -= 2;
    }
    if(!(x & 0x80000000u))
    {
      x <<= 1;
      r -= 1;
    }
    return r;
  }

}  // namespace utility

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_UTILITY_FUNCTIONS_HH
