//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ProcTblElement.cc,v 1.6 2001-07-11 10:08:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
























