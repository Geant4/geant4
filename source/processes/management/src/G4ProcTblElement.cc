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
// G4ProcTblElement class implementation
//
// Author: H.Kurashige, 4 August 1998
// --------------------------------------------------------------------

#include "G4ProcTblElement.hh"

// --------------------------------------------------------------------
// Default constructor
//
G4ProcTblElement::G4ProcTblElement()
{
}

// --------------------------------------------------------------------
G4ProcTblElement::G4ProcTblElement(G4VProcess* aProcess)
  : pProcess(aProcess)
{
  pProcMgrVector = new G4ProcMgrVector();
}

// --------------------------------------------------------------------
// Copy constructor
//
G4ProcTblElement::G4ProcTblElement(const G4ProcTblElement& right)
{
  *this = right;
}

// --------------------------------------------------------------------
// Destructor
//
G4ProcTblElement::~G4ProcTblElement()
{
  if (pProcMgrVector != nullptr)
  {
    pProcMgrVector->clear();
    delete pProcMgrVector;
  }
}

// --------------------------------------------------------------------
G4ProcTblElement& G4ProcTblElement::operator=(const G4ProcTblElement& right)
{
  if (this != &right)
  {
    pProcess = right.pProcess;

    // copy all contents in  pProcMgrVector
    if (pProcMgrVector != nullptr)
    {
      pProcMgrVector->clear();
      delete pProcMgrVector;
    }
    pProcMgrVector = new G4ProcMgrVector();
    for (auto i = right.pProcMgrVector->cbegin();
              i!= right.pProcMgrVector->cend(); ++i)
    {
      pProcMgrVector->push_back(*i);
    }
  }
  return *this;
}

// --------------------------------------------------------------------
G4bool G4ProcTblElement::operator==(const G4ProcTblElement& right) const
{
  return (this == &right);
}

// --------------------------------------------------------------------
G4bool G4ProcTblElement::operator!=(const G4ProcTblElement& right) const
{
  return (this != &right);
}
