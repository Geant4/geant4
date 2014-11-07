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
// $Id: $
//
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//    A singleton to maintain the list of G4VBiasingOperation objects.
//
//      ----------------G4BiasingOperationManager ----------------
//
// Author: M.Verderi (LLR), November 2013
//
// --------------------------------------------------------------------

#ifndef G4BiasingOperationManager_hh
#define G4BiasingOperationManager_hh 1

class G4VBiasingOperation;
#include <map>
#include <vector>
#include "G4Cache.hh"
#include "G4ThreadLocalSingleton.hh"

//This class is a thread-local singleton 
class G4BiasingOperationManager {
    friend class G4ThreadLocalSingleton<G4BiasingOperationManager>;
public:
  static G4BiasingOperationManager*                  GetInstance();
  const std::vector< G4VBiasingOperation* > GetBiasingOperations() {return fBiasingOperationVector.Get();}
  G4VBiasingOperation*                       GetBiasingOperation(std::size_t optionID);

public:
  ~G4BiasingOperationManager();
  std::size_t Register(G4VBiasingOperation*);
  
  
private:
  G4BiasingOperationManager();
  static G4VectorCache<G4VBiasingOperation*> fBiasingOperationVector;
  static G4MapCache<G4VBiasingOperation*,std::size_t > fBiasingOperationIDtoPointerMap;
};

#endif
