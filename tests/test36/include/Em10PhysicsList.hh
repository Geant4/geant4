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
// $Id: Em10PhysicsList.hh,v 1.4 2007-06-21 15:06:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Em10PhysicsList_h
#define Em10PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4PhotoElectricEffect;
class G4ComptonScattering;
class G4GammaConversion;

class G4eMultipleScattering;

class G4PAIonisation ;
class G4ForwardXrayTR ;
class G4eIonisation;
class G4eBremsstrahlung;
class G4eplusAnnihilation;

class G4MuIonisation;
class G4MuBremsstrahlung;
class G4MuPairProduction;

class G4hIonisation;

class Em10StepCut;

class Em10DetectorConstruction;
// class ALICEDetectorConstruction;
class Em10PhysicsListMessenger;
class G4ProductionCuts;


class Em10PhysicsList: public G4VModularPhysicsList  // G4VUserPhysicsList
{
  public:
    Em10PhysicsList( Em10DetectorConstruction*);
  // Em10PhysicsList( ALICEDetectorConstruction*);
   ~Em10PhysicsList();

  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
 
    void SetCuts();

  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBarions();

  protected:

    // these methods Construct physics processes and register them

    void AddParameterisation();
    void ConstructGeneral();
    void ConstructEM();

  public:

    void SetGammaCut(G4double);
    void SetElectronCut(G4double);

    void SetRegGammaCut(G4double    cut){fGammaCut    = cut;};
    void SetRegElectronCut(G4double cut){fElectronCut = cut;};
    void SetRegPositronCut(G4double cut){fPositronCut = cut;};

    void SetRadiatorCuts();
    void SetDetectorCuts();

    void SetMaxStep(G4double);
    void SetMinElectronEnergy(G4double E){fMinElectronEnergy=E;};     
    void SetMinGammaEnergy(G4double E)   {fMinGammaEnergy=E;};       
  void SetXTRModel(G4String m)   {fXTRModel = m; G4cout<<fXTRModel<<G4endl;};       

  public:   

    G4double MaxChargedStep;

  private:

    G4PhotoElectricEffect* thePhotoElectricEffect;
    G4ComptonScattering*   theComptonScattering;
    G4GammaConversion*     theGammaConversion;

    G4eMultipleScattering*  theeminusMultipleScattering;
    G4eIonisation*         theeminusIonisation;
    G4eBremsstrahlung*     theeminusBremsstrahlung;

    G4ForwardXrayTR*       fForwardXrayTR ;

    G4eMultipleScattering*  theeplusMultipleScattering;
    G4eIonisation*         theeplusIonisation;
    G4eBremsstrahlung*     theeplusBremsstrahlung;
    G4eplusAnnihilation*   theeplusAnnihilation;

    Em10StepCut* theeminusStepCut ;
    Em10StepCut* theeplusStepCut ;

    G4double cutForGamma;
    G4double cutForElectron, cutForPositron;

    Em10DetectorConstruction* pDet;
  //  ALICEDetectorConstruction* apDet;
    Em10PhysicsListMessenger* physicsListMessenger;

    G4double fMinElectronEnergy;      // minimalEnergy of produced electrons
    G4double fMinGammaEnergy; 
        // minimalEnergy of scattered photons
    G4ProductionCuts* fRadiatorCuts;
    G4ProductionCuts* fDetectorCuts;
    G4double fElectronCut, fGammaCut, fPositronCut;
    G4String fXTRModel;
};

#endif



