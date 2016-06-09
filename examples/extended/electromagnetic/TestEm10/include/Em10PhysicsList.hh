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
// $Id: Em10PhysicsList.hh,v 1.4 2003/08/28 09:36:12 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
//

#ifndef Em10PhysicsList_h
#define Em10PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4PhotoElectricEffect;
class G4ComptonScattering;
class G4GammaConversion;

class G4MultipleScattering52;

class G4PAIonisation ;
class G4ForwardXrayTR ;
class G4eIonisation52;
class G4eBremsstrahlung52;
class G4eplusAnnihilation;

class G4MuIonisation52;
class G4MuBremsstrahlung52;
class G4MuPairProduction52;

class G4hIonisation52;

class Em10StepCut;

class Em10DetectorConstruction;
class Em10PhysicsListMessenger;


class Em10PhysicsList: public G4VUserPhysicsList
{
  public:
    Em10PhysicsList( Em10DetectorConstruction*);
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

    void SetMaxStep(G4double);

  public:   

    G4double MaxChargedStep;

  private:

    G4PhotoElectricEffect* thePhotoElectricEffect;
    G4ComptonScattering*   theComptonScattering;
    G4GammaConversion*     theGammaConversion;

    G4MultipleScattering52*  theeminusMultipleScattering;
    G4eIonisation52*         theeminusIonisation;
    G4eBremsstrahlung52*     theeminusBremsstrahlung;

    G4PAIonisation*        fPAIonisation ;
    G4ForwardXrayTR*       fForwardXrayTR ;

    G4MultipleScattering52*  theeplusMultipleScattering;
    G4eIonisation52*         theeplusIonisation;
    G4eBremsstrahlung52*     theeplusBremsstrahlung;
    G4eplusAnnihilation*   theeplusAnnihilation;

    Em10StepCut* theeminusStepCut ;
    Em10StepCut* theeplusStepCut ;

    G4double cutForGamma;
    G4double cutForElectron;

    Em10DetectorConstruction* pDet;
    Em10PhysicsListMessenger* physicsListMessenger;
};

#endif



