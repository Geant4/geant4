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
// 26 August 2023   V.Ivanchenko
//
// This is a singleton class to store shared G4ElementData
//

#ifndef G4ElementDataRegistry_h
#define G4ElementDataRegistry_h 1

#include <vector>
#include "G4ElementData.hh"

class G4ElementDataRegistry
{
 public:

  static G4ElementDataRegistry* Instance();

  ~G4ElementDataRegistry();

  void RegisterMe(G4ElementData* p);

  void RemoveMe(G4ElementData* p);

  const std::vector<G4ElementData*>& GetElementData() const {
    return elmdata;
  }

  G4ElementData* GetElementDataByName(const G4String&);

 private:

  G4ElementDataRegistry();

  static G4ElementDataRegistry* instance;

  std::vector<G4ElementData*> elmdata;
};

#endif
