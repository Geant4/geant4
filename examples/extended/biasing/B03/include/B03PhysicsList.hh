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
/// \file biasing/B03/include/B03PhysicsList.hh
/// \brief Definition of the B03PhysicsList class
//
//
#ifndef B03PhysicsList_h
#define B03PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4GeometrySampler.hh"

#include <vector>

// taken from Tst12PhysicsList

class B03PhysicsList: public G4VUserPhysicsList
{
  public:
    B03PhysicsList(G4String);
    virtual ~B03PhysicsList();

  public: 
    void AddParallelWorldName(G4String& pname)
         {fParaWorldName.push_back(pname);}

    // void AddBiasing(G4GeometrySampler *mgs, G4String& pname) 
    //      {fGeomSampler = mgs; fBiasWorldName = pname;}

  protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // 
    virtual void SetCuts();
    
  protected:
  // these methods Construct physics processes and register them
    virtual void ConstructGeneral();
    virtual void ConstructEM();
    virtual void ConstructHad();
    virtual void ConstructLeptHad();

    void AddScoringProcess(); 
    void AddBiasingProcess(); 

 //
    void  ConstructAllBosons();
    void  ConstructAllLeptons();
    void  ConstructAllMesons();
    void  ConstructAllBaryons();
    void  ConstructAllIons();
    void  ConstructAllShortLiveds();

  private:
    std::vector<G4String>  fParaWorldName; 
  G4String fBiasWorldName;
  // G4GeometrySampler* fGeomSampler;

};

#endif

