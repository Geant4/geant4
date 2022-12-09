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
#include "UtilityFunctions.hh"

#include <memory>
#include <cstdio>

#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace utility
{
  std::vector<G4String>& Split(const G4String& str, char delim,
                               std::vector<G4String>& elems)
  {
    std::stringstream ss(str);
    std::string item;
    while(std::getline(ss, item, delim))
    {
      if(!item.empty())
        elems.push_back((G4String) item);
    }
    return elems;
  }

  std::vector<G4String> Split(const G4String& str, char delim)
  {
    std::vector<G4String> elems;
    Split(str, delim, elems);
    return elems;
  }

  // return element ii
  G4String Get_seperated_element(const G4String& str, char delim, G4int ii)
  {
    G4int delims = 0;
    G4String ss;
    for(char c : str)
    {
      if(delims > ii) {
        return ss;}
      if(c == delim) {
        delims++;
      } else if(delims == ii) {
        ss += c;}
    }
    return ss;
  }

  // return first four strings
  std::array<G4String, 4> Get_four_elements(const G4String& str, char delim)
  {
    G4int delims = 0;
    std::array<G4String, 4> arr;
    for(char c : str)
    {
      if(delims >= 4)
        return arr;
      if(c == delim)
        delims++;
      else
        arr.at(delims) += c;
    }
    return arr;
  }

  G4bool Path_exists(const G4String& fname)
  {
    if(FILE* file = std::fopen(fname, "r"))
    {
      fclose(file);
      return true;
    }
    else
    {
      return false;
    }
  }

  // Memory safe multi-platform getcwd
  // http://stackoverflow.com/questions/2869594/
  G4String Getcwd()
  {
    const size_t chunkSize = 255;
    const int maxChunks    = 10240;  // 2550 KiBs of current path is plenty

    char stackBuffer[chunkSize];  // Stack buffer for the "normal" case
    if(::getcwd(stackBuffer, sizeof(stackBuffer)) != nullptr)
      return stackBuffer;
    if(errno != ERANGE)
    {
      // It's not ERANGE, so we don't know how to handle it
      throw std::runtime_error("Cannot determine the current path.");
      // Of course you may choose a different error reporting method
    }
    // Ok, the stack buffer isn't long enough; fallback to heap allocation
    for(int chunks = 2; chunks < maxChunks; chunks++)
    {
      // With boost use scoped_ptr; in C++0x, use unique_ptr
      // If you want to be less C++ but more efficient
      // you may want to use realloc
      std::unique_ptr<char[]> cwd(new char[chunkSize * chunks]);
      if(::getcwd(cwd.get(), chunkSize * chunks) != nullptr)
        return cwd.get();
      if(errno != ERANGE)
      {
        // It's not ERANGE, so we don't know how to handle it
        throw std::runtime_error("Cannot determine the current path.");
        // Of course you may choose a different error reporting method
      }
    }
    throw std::runtime_error("Cannot determine the current path; the path is "
                             "apparently unreasonably long");
  }

  G4double Min(const G4ThreeVector& v)
  {
    return std::min(std::min(v.x(), v.y()), v.z());
  }
}  // namespace utility

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
