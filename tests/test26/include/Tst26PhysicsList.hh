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
// $Id: Tst26PhysicsList.hh,v 1.6 2003-03-04 19:08:22 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
// 14-02-03 Make G4ProductionCuts to be members of the class (V.Ivanchenko)
// 19-02-03 Rename G4ProductionCuts (V.Ivanchenko)
// 04-03-03 Define default EM module (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//

//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst26PhysicsList_h
#define Tst26PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class Tst26PhysicsListMessenger;
class G4ProductionCuts;
class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst26PhysicsList: public G4VModularPhysicsList
{
  public:
    Tst26PhysicsList();
   ~Tst26PhysicsList();


    virtual void ConstructParticle();
    virtual void ConstructProcess();
    void AddPhysicsList(const G4String& name);

    void SetCuts();
    void SetCutForWorld(G4double);
    void SetCutForVertexDetector(G4double);
    void SetCutForMuonDetector(G4double);

  private:
    G4double cutForWorld;
    G4double cutForVertexDetector;
    G4double cutForMuonDetector;

//    typedef G4std::vector<G4VPhysicsConstructor*> G4PhysConstVector;
//    G4PhysConstVector* physicsVector;

    G4VPhysicsConstructor*  emPhysicsList;
    G4VPhysicsConstructor*  generalPhysicsList;
    G4VPhysicsConstructor*  particleList;
    G4String emName;

    Tst26PhysicsListMessenger* pMessenger;
    G4ProductionCuts* vertexDetectorCuts;
    G4ProductionCuts* muonDetectorCuts;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

