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
// $Id: G4ProcTblElement.cc 71231 2013-06-12 13:06:28Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	4th Aug 1998, H.Kurashige
// ------------------------------------------------------------
//   Use STL vector instead of RW vector    1. Mar 00 H.Kurashige
//
#include "G4ProcTblElement.hh"


// default constructor ////////////////////////
G4ProcTblElement::G4ProcTblElement()
   :pProcess(0),pProcMgrVector(0)
{
}

//////////////////////////
G4ProcTblElement::G4ProcTblElement(G4VProcess* aProcess):
  pProcess(aProcess)
{
  pProcMgrVector = new G4ProcMgrVector();
}

// copy constructor //////////////////////////
G4ProcTblElement::G4ProcTblElement(const G4ProcTblElement &right)
  :pProcess(0),pProcMgrVector(0)
{
  *this = right;
}


// destructor ////////////////////////
G4ProcTblElement::~G4ProcTblElement()
{
  if (pProcMgrVector != 0) {
    pProcMgrVector->clear();
    delete pProcMgrVector;
  }
}


//////////////////////////
G4ProcTblElement & G4ProcTblElement::operator=(const G4ProcTblElement &right)
{
  if (this != &right) {
    pProcess       = right.pProcess;
    // copy all contents in  pProcMgrVector
    if (pProcMgrVector != 0) {
      pProcMgrVector->clear();
      delete pProcMgrVector;
    }
    pProcMgrVector = new G4ProcMgrVector();
    G4ProcMgrVector::iterator i;
    for (i = right.pProcMgrVector->begin(); i!= right.pProcMgrVector->end(); ++i) {
      pProcMgrVector->push_back(*i);
    }
  }
  return *this;
}


//////////////////////////
G4int G4ProcTblElement::operator==(const G4ProcTblElement &right) const
{
  return (this == &right);
}


//////////////////////////
G4int G4ProcTblElement::operator!=(const G4ProcTblElement &right) const
{
  return (this != &right);
}
























