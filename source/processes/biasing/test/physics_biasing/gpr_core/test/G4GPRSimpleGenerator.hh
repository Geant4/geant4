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
// $Id: G4GPRSimpleGenerator.hh,v 1.1 2007-08-02 18:12:06 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPRSIMPLeGENERATOR_HH
#define G4GPRSIMPLEGENERATOR_HH

#include "G4GPRElementSuperStore.hh"
#include "G4GPRUtils.hh"

class G4GPRSimpleGenerator {

public:

  G4GPRSimpleGenerator()
  {
    fSuperStore = G4GPRElementSuperStore::Instance();
    //    pCacheManager = G4GPRCacheManagerSuperStore::Instance();
  }

  template <typename T, typename Result>
  void Generate(Result*& result) 
  {
    result = new Result;

    G4GPRElementStoreT<T>* store = fSuperStore;
    
    G4GPRUtils::Operator(result, store);
  }

private:

  G4GPRElementSuperStore* fSuperStore;

};

#endif
