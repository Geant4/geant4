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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRPhysicsList.hh
//   Header file of the Gorad Physics List
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRPhysicsList_H
#define GRPhysicsList_H 1

#include "G4VModularPhysicsList.hh"
class G4PhysListFactory;
class GRPhysicsListMessenger;
class G4Region;
class G4ProductionCuts;
 
#include <map>

class GRPhysicsList : public G4VModularPhysicsList
{
  public:
    GRPhysicsList();
    virtual ~GRPhysicsList();
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    virtual void SetCuts();

  private:
    G4String PLName;
    G4VModularPhysicsList* physList;
    G4PhysListFactory* factory;
    GRPhysicsListMessenger* messenger;

  public:
    const G4String& GetPLName()
    { return PLName; }

  private:
    G4String EM_opt;   // EM physics option
    G4String Had_opt;  // Hadronic physics option
    G4bool addHP;      // add Neutron_HP
    G4bool addRDM;     // add Radioactive Decay Module
    G4bool addRMC;      // add Reverse Monte Calro
    G4bool addOptical; // add optical physics
    G4int stepLimit_opt; // Step limiter option (0:charged, 1:neutral, 2:all, 3:e+/e-)
    std::map<G4Region*,G4double> localStepLimits; // map of region name and limit value
    G4double globalCuts[4];  // for e-, e+ gamma, proton
    std::map<G4Region*,G4ProductionCuts*> localCuts; // map of region name and cuts

  public:
    void SetEM(G4String& val) { EM_opt = val; }
    const G4String& GetEM() const { return EM_opt; }
    void SetHad(G4String& val) { Had_opt = val; }
    const G4String& GetHad() const { return Had_opt; }
    void AddHP(G4bool val = true) { addHP = val; }
    G4bool IfHP() const { return addHP; }
    void AddRDM(G4bool val = true) { addRDM = val; }
    G4bool IfRDM() const { return addRDM; }
    void AddRMC(G4bool val = true) { addRMC = val; }
    G4bool IfRMC() const { return addRMC; }
    void AddOptical(G4bool val = true) { addOptical = val; }
    G4bool IfOptical() const { return addOptical; }
    void AddStepLimit(G4int val = 0) { stepLimit_opt = val; }
    G4int IfStepLimit() const { return stepLimit_opt; }
    void SetGlobalStepLimit(G4double);
    G4double GetGlobalStepLimit() const;
    G4Region* SetLocalStepLimit(const G4String&,G4double);
    G4double GetLocalStepLimit(const G4String&) const;
    void SetGlobalCuts(G4double);
    G4double GetGlobalCuts() const { return GetGlobalCut(0); }
    void SetGlobalCut(G4int, G4double); 
    G4double GetGlobalCut(G4int i) const { return globalCuts[i]; }
    G4Region* SetLocalCuts(const G4String& reg,G4double val)
    {
      G4Region* regPtr = nullptr;
      for(G4int i=0; i<4; i++)
      {
        regPtr = SetLocalCut(reg,i,val);
        if(!regPtr) return regPtr;
      }
      return regPtr;
    }
    G4double GetLocalCuts(const G4String& reg) const { return GetLocalCut(reg,0); }
    G4Region* SetLocalCut(const G4String&,G4int,G4double);
    G4double GetLocalCut(const G4String&,G4int) const;

  private:
    void GeneratePLName();
    void GeneratePL();
    G4Region* FindRegion(const G4String&) const;

  private:
    G4bool applyGeomImpBias = false;

  public:
    void ApplyGeomImpBias(G4bool val = true)
    { applyGeomImpBias = val; }
};

#endif

