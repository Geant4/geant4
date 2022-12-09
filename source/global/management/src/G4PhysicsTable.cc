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
// G4PhysicsTable class implementation
//
// Author: G.Cosmo, 2 December 1995
//         First implementation based on object model
// Revisions:
// - 1st March 1996, K.Amako: modified
// - 24th February 2001, H.Kurashige: migration to STL vectors
// --------------------------------------------------------------------

#include <fstream>
#include <iomanip>
#include <iostream>

#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsVectorType.hh"

// --------------------------------------------------------------------
G4PhysicsTable::G4PhysicsTable(std::size_t cap)
{
  reserve(cap);
  vecFlag.reserve(cap);
}

// --------------------------------------------------------------------
G4PhysicsTable::~G4PhysicsTable()
{
  G4PhysCollection::clear();
  vecFlag.clear();
}

// --------------------------------------------------------------------
void G4PhysicsTable::resize(std::size_t siz, G4PhysicsVector* vec)
{
  G4PhysCollection::resize(siz, vec);
  vecFlag.resize(siz, true);
}

// --------------------------------------------------------------------
G4bool G4PhysicsTable::StorePhysicsTable(const G4String& fileName, G4bool ascii)
{
  std::ofstream fOut;

  // open output file
  if(!ascii)
  {
    fOut.open(fileName, std::ios::out | std::ios::binary);
  }
  else
  {
    fOut.open(fileName, std::ios::out);
  }

  // check if the file has been opened successfully
  if(!fOut.is_open())
  {
#ifdef G4VERBOSE
    G4cerr << "G4PhysicsTable::StorePhysicsTable():";
    G4cerr << " Cannot open file: " << fileName << G4endl;
#endif
    fOut.close();
    return false;
  }

  // Number of elements
  std::size_t tableSize = size();
  if(!ascii)
  {
    fOut.write((char*) (&tableSize), sizeof tableSize);
  }
  else
  {
    fOut << tableSize << G4endl;
  }

  // Physics Vector
  for(const auto itr : *this)
  {
    G4int vType = itr->GetType();
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
G4bool G4PhysicsTable::ExistPhysicsTable(const G4String& fileName) const
{
  std::ifstream fIn;
  G4bool value = true;
  // open input file
  fIn.open(fileName, std::ios::in);

  // check if the file has been opened successfully
  if(!fIn)
  {
    value = false;
  }
  fIn.close();
  return value;
}

// --------------------------------------------------------------------
G4bool G4PhysicsTable::RetrievePhysicsTable(const G4String& fileName,
                                            G4bool ascii, G4bool spline)
{
  std::ifstream fIn;
  // open input file
  if(ascii)
  {
    fIn.open(fileName, std::ios::in | std::ios::binary);
  }
  else
  {
    fIn.open(fileName, std::ios::in);
  }

  // check if the file has been opened successfully
  if(!fIn.is_open())
  {
#ifdef G4VERBOSE
    G4cerr << "G4PhysicsTable::RetrievePhysicsTable():";
    G4cerr << " Cannot open file: " << fileName << G4endl;
#endif
    fIn.close();
    return false;
  }

  // clear
  clearAndDestroy();

  // Number of elements
  std::size_t tableSize = 0;
  if(!ascii)
  {
    fIn.read((char*) (&tableSize), sizeof tableSize);
  }
  else
  {
    fIn >> tableSize;
  }
  reserve(tableSize);
  vecFlag.clear();

  // Physics Vector
  for(std::size_t idx = 0; idx < tableSize; ++idx)
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
    G4PhysicsVector* pVec = CreatePhysicsVector(vType, spline);
    if(pVec == nullptr)
    {
#ifdef G4VERBOSE
      G4cerr << "G4PhysicsTable::RetrievePhysicsTable():";
      G4cerr << " Illegal Physics Vector type: " << vType << " in: ";
      G4cerr << fileName << G4endl;
#endif
      fIn.close();
      return false;
    }

    if(!(pVec->Retrieve(fIn, ascii)))
    {
#ifdef G4VERBOSE
      G4cerr << "G4PhysicsTable::RetrievePhysicsTable():";
      G4cerr << " Rrror in retreiving " << idx
             << "-th Physics Vector from file: ";
      G4cerr << fileName << G4endl;
#endif
      fIn.close();
      return false;
    }

    // add a PhysicsVector to this PhysicsTable
    G4PhysCollection::push_back(pVec);
    vecFlag.push_back(true);
  }
  fIn.close();
  return true;
}

// --------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, G4PhysicsTable& right)
{
  // Printout Physics Vector
  std::size_t i = 0;
  for(auto itr = right.cbegin(); itr != right.cend(); ++itr)
  {
    out << std::setw(8) << i << "-th Vector   ";
    out << ": Type    " << G4int((*itr)->GetType());
    out << ": Flag    ";
    if(right.GetFlag(i))
    {
      out << " T";
    }
    else
    {
      out << " F";
    }
    out << G4endl;
    out << *(*itr);
    ++i;
  }
  out << G4endl;
  return out;
}

// --------------------------------------------------------------------
void G4PhysicsTable::ResetFlagArray()
{
  size_t tableSize = G4PhysCollection::size();
  vecFlag.clear();
  for(std::size_t idx = 0; idx < tableSize; ++idx)
  {
    vecFlag.push_back(true);
  }
}

// --------------------------------------------------------------------
G4PhysicsVector* G4PhysicsTable::CreatePhysicsVector(G4int type, G4bool spline)
{
  G4PhysicsVector* pVector = nullptr;
  switch(type)
  {
    case T_G4PhysicsLinearVector:
      pVector = new G4PhysicsLinearVector(spline);
      break;

    case T_G4PhysicsLogVector:
      pVector = new G4PhysicsLogVector(spline);
      break;

    default:
      pVector = new G4PhysicsVector(spline);
      break;
  }
  return pVector;
}
