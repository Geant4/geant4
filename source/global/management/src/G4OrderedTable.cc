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
// G4OrderedTable class implementation
//
// Author: M.Maire (LAPP), September 1996
// Revisions: H.Kurashige (Kobe Univ.), January-September 2001
// --------------------------------------------------------------------

#include "G4OrderedTable.hh"
#include "G4DataVector.hh"

#include <fstream>
#include <iomanip>
#include <iostream>

// --------------------------------------------------------------------
G4OrderedTable::G4OrderedTable(std::size_t cap)
  : std::vector<G4DataVector*>(cap, (G4DataVector*) nullptr)
{}

// --------------------------------------------------------------------
void G4OrderedTable::clearAndDestroy()
{
  G4DataVector* a = nullptr;
  while(!empty())
  {
    a = back();
    pop_back();
    for(auto i = cbegin(); i != cend(); ++i)
    {
      if(*i == a)
      {
        erase(i);
        --i;
      }
    }

    delete a;
  }
}

// --------------------------------------------------------------------
G4bool G4OrderedTable::Store(const G4String& fileName, G4bool ascii)
{
  std::ofstream fOut;

  // open output file //
  if(!ascii)
  {
    fOut.open(fileName, std::ios::out | std::ios::binary);
  }
  else
  {
    fOut.open(fileName, std::ios::out);
  }

  // check if the file has been opened successfully
  if(!fOut)
  {
#ifdef G4VERBOSE
    G4cerr << "G4OrderedTable::::Store():";
    G4cerr << " Cannot open file: " << fileName << G4endl;
#endif
    fOut.close();
    return false;
  }

  G4int tableSize = G4int(size());  // Number of elements
  if(!ascii)
  {
    fOut.write((char*) (&tableSize), sizeof tableSize);
  }
  else
  {
    fOut << tableSize << G4endl;
  }

  G4int vType = G4DataVector::T_G4DataVector;  // Data Vector
  for(const auto itr : *this)
  {
    if(!ascii)
    {
      fOut.write((char*) (&vType), sizeof vType);
    }
    else
    {
      fOut << vType << G4endl;
    }
    itr->Store(fOut, ascii);
  }
  fOut.close();
  return true;
}

// --------------------------------------------------------------------
G4bool G4OrderedTable::Retrieve(const G4String& fileName, G4bool ascii)
{
  std::ifstream fIn;
  // open input file //
  if(!ascii)
  {
    fIn.open(fileName, std::ios::in | std::ios::binary);
  }
  else
  {
    fIn.open(fileName, std::ios::in);
  }

  // check if the file has been opened successfully
  if(!fIn)
  {
#ifdef G4VERBOSE
    G4cerr << "G4OrderedTable::Retrieve():";
    G4cerr << " Cannot open file: " << fileName << G4endl;
#endif
    fIn.close();
    return false;
  }

  // clear
  clearAndDestroy();

  // Number of elements
  G4int tableSize = 0;
  if(!ascii)
  {
    fIn.read((char*) (&tableSize), sizeof tableSize);
  }
  else
  {
    fIn >> tableSize;
  }
  if(tableSize <= 0)
  {
#ifdef G4VERBOSE
    G4cerr << "G4OrderedTable::Retrieve():";
    G4cerr << " Invalid table size: " << tableSize << G4endl;
#endif
    return false;
  }
  reserve(tableSize);

  // Physics Vector
  for(G4int idx = 0; idx < tableSize; ++idx)
  {
    G4int vType = 0;
    if(!ascii)
    {
      fIn.read((char*) (&vType), sizeof vType);
    }
    else
    {
      fIn >> vType;
    }
    if(vType != G4DataVector::T_G4DataVector)
    {
#ifdef G4VERBOSE
      G4cerr << "G4OrderedTable::Retrieve():";
      G4cerr << " Illegal Data Vector type: " << vType << " in  ";
      G4cerr << fileName << G4endl;
#endif
      fIn.close();
      return false;
    }

    auto* pVec = new G4DataVector;

    if(!(pVec->Retrieve(fIn, ascii)))
    {
#ifdef G4VERBOSE
      G4cerr << "G4OrderedTable::Retrieve(): ";
      G4cerr << " Error in retreiving " << idx
             << "-th Physics Vector from file: ";
      G4cerr << fileName << G4endl;
#endif
      fIn.close();
      delete pVec;
      return false;
    }

    // add a PhysicsVector to this OrderedTable
    push_back(pVec);
  }
  fIn.close();
  return true;
}

// --------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, G4OrderedTable& right)
{
  // Printout Data Vector
  std::size_t i = 0;
  for(const auto itr : right)
  {
    out << std::setw(8) << i << "-th Vector   ";
    out << ": Type    " << G4DataVector::T_G4DataVector << G4endl;
    out << *itr;
    i += 1;
  }
  out << G4endl;
  return out;
}
