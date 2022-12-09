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
//  G4Hashing.cc
//  Geant4
//
//  Created by Mathieu Karamitros
//
//

#include "DNAHashing.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace G4::hashing
  {
    namespace crc32
    {
      uint32_t Hash(const char* str, size_t len)
      {
        uint32_t remainder = 0xFFFFFFFF;
        for(size_t idx = 0; idx < len; ++idx)
        {
          remainder =
            (remainder >> 8) ^ fCrc_table[(remainder ^ str[idx]) & 0x000000FF];
        }
        return remainder ^ 0xFFFFFFFF;
      }

      uint32_t Hash(const std::string& str)
      {
        return Hash(str.c_str(), str.size());
      }
    }  // namespace crc32

    namespace fnv
    {
      size_t Hash(const std::string& str)
      {
        size_t hash = fnv_offset_basis;
        for(size_t i = 0; i < str.length(); i++)
        {
          hash = hash ^ (str[i]);   // xor  the low 8 bits
          hash = hash * fnv_prime;  // multiply by the magic number
        }
        return hash;
      }
    }  // namespace fnv

    namespace larson
    {
      size_t Hash(const char* str, unsigned int seed)
      {
        size_t hash = seed;
        while(*str)
        {
          hash = hash * 101 + *str++;
        }
        return hash;
      }

      size_t Hash(std::string str, unsigned int seed)
      {
        return Hash(str.c_str(), seed);
      }
    }  // namespace larson
  }  // namespace G4

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
