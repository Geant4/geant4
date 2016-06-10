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
// $Id: G4ProcessVector.cc 71231 2013-06-12 13:06:28Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//  This class is a container for pointers to physics process objects.
// ------------------------------------------------------------

#include "G4VProcess.hh"
#include "G4ProcessVector.hh"

/////// Constructors
//
G4ProcessVector::G4ProcessVector()
{
  pProcVector = new G4ProcVector();
}

G4ProcessVector::G4ProcessVector(size_t siz)
{
  pProcVector = new G4ProcVector(siz);
}

G4ProcessVector::G4ProcessVector(const G4ProcessVector& right)
  :pProcVector(0)
{
  pProcVector = new G4ProcVector();
  
  // copy all contents in  pProcVector
  //
  G4ProcVector::iterator i;
  for (i = right.pProcVector->begin(); i!= right.pProcVector->end(); ++i){
    pProcVector->push_back(*i);
  }
}

/////// destructor
//
G4ProcessVector::~G4ProcessVector()
{
  // delete pProcVector
  //
  if (pProcVector != 0 ){
    pProcVector->clear();
    delete pProcVector;
  }
}

////// assignment oeprator
//
G4ProcessVector & G4ProcessVector::operator=(const G4ProcessVector & right)
{
  if (this != &right){
    // delete pProcVector
    if (pProcVector != 0 ){
      pProcVector->clear();
      delete pProcVector;
    }
    
    pProcVector = new G4ProcVector();
    
    // copy all contents in  pProcVector
    //
    G4ProcVector::iterator i;
    for (i = right.pProcVector->begin(); i!= right.pProcVector->end(); ++i){
      pProcVector->push_back(*i);
    }
  }
  return *this;
}
