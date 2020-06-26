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
// G4ProcessVector class implementation
//
// Authors: G.Cosmo, H.Kurashige - 1998
// --------------------------------------------------------------------

#include "G4VProcess.hh"
#include "G4ProcessVector.hh"

// --------------------------------------------------------------------
/////// Constructors
//
G4ProcessVector::G4ProcessVector()
{
  pProcVector = new G4ProcVector();
}

G4ProcessVector::G4ProcessVector(std::size_t siz)
{
  pProcVector = new G4ProcVector(siz);
}

G4ProcessVector::G4ProcessVector(const G4ProcessVector& right)
  : pProcVector(nullptr)
{
  pProcVector = new G4ProcVector();
  
  // copy all contents in  pProcVector
  //
  for (auto i=right.pProcVector->cbegin(); i!=right.pProcVector->cend(); ++i)
  {
    pProcVector->push_back(*i);
  }
}

// --------------------------------------------------------------------
/////// destructor
//
G4ProcessVector::~G4ProcessVector()
{
  // delete pProcVector
  //
  if (pProcVector != nullptr )
  {
    pProcVector->clear();
    delete pProcVector;
  }
}

// --------------------------------------------------------------------
////// assignment operator
//
G4ProcessVector& G4ProcessVector::operator=(const G4ProcessVector& right)
{
  if (this != &right)
  {
    // delete pProcVector
    if (pProcVector != nullptr )
    {
      pProcVector->clear();
      delete pProcVector;
    }
    
    pProcVector = new G4ProcVector();
    
    // copy all contents in  pProcVector
    //
    for (auto i=right.pProcVector->cbegin(); i!=right.pProcVector->cend(); ++i)
    {
      pProcVector->push_back(*i);
    }
  }
  return *this;
}

// --------------------------------------------------------------------
//
std::size_t G4ProcessVector::index(G4VProcess* aProcess) const
{
  std::size_t j=0;
  for (auto it=pProcVector->cbegin();it!=pProcVector->cend(); ++j, ++it)
  {
    if (*(*it)==*aProcess) return j;
  }
  return -1;
}
 
// --------------------------------------------------------------------
//
G4bool G4ProcessVector::contains(G4VProcess* aProcess) const
{
  for (auto it=pProcVector->cbegin();it!=pProcVector->cend(); ++it)
  {
    if (*(*it)==*aProcess) return true;
  }
  return false;
}

// --------------------------------------------------------------------
//
G4bool G4ProcessVector::insertAt(G4int i, G4VProcess* aProcess)
{
  if ( (i<0) || (i>G4int(pProcVector->size())) ) return false;
  if (i==G4int(pProcVector->size()))
  {
    pProcVector->push_back(aProcess);
  }
  else
  {
    auto it=pProcVector->cbegin();
    for(G4int j=0; j!=i; ++j) ++it;
    pProcVector->insert(it, aProcess);
  }
  return true;
}

// --------------------------------------------------------------------
//
G4VProcess* G4ProcessVector::removeAt(G4int i)
{
  auto it=pProcVector->cbegin();
  for(std::size_t j=0; j<pProcVector->size() && G4int(j)<i; ++j) ++it;
  G4VProcess* rValue = *it;
  pProcVector->erase(it);
  return rValue;
}
