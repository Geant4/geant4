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
// $Id: $
//
// File:    G4GDecay3.hh
// Author:  Dennis Wright (SLAC)
// Date:    19 April 2013
//
// Description: three-body phase space momentum generator based on 
//              GDECA3 of Geant3
//
#ifndef G4GDecay3_h
#define G4GDecay3_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>


class G4GDecay3 {

  public:
    G4GDecay3() {;}
    G4GDecay3(const G4double& pMass, const G4double& dMass0,
              const G4double& dMass1, const G4double& dMass2);
    ~G4GDecay3() {;}
  
    std::vector<G4ThreeVector> GetThreeBodyMomenta();

  private:
    G4bool CalculateMomentumMagnitudes();

    G4int loopMax;

    G4double parentMass;
    G4double mDaughter0;
    G4double mDaughter1;
    G4double mDaughter2;

    G4double pDaughter0;
    G4double pDaughter1;
    G4double pDaughter2;
};        

#endif

