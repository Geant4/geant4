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
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSOpticalPhysics_h
#define WLSOpticalPhysics_h 1

#include "globals.hh"

#include "G4OpWLS.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"

#include "G4OpMieHG.hh"
#include "G4OpRayleigh.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4VPhysicsConstructor.hh"

class WLSOpticalPhysics : public G4VPhysicsConstructor
{
  public:

    WLSOpticalPhysics(G4bool toggle=true);
    virtual ~WLSOpticalPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

    G4OpWLS* GetWLSProcess() {return theWLSProcess;}
    G4Cerenkov* GetCerenkovProcess() {return theCerenkovProcess;}
    G4Scintillation* GetScintillationProcess() {return theScintProcess;}
    G4OpAbsorption* GetAbsorptionProcess() {return theAbsorptionProcess;}
    G4OpRayleigh* GetRayleighScatteringProcess() {return theRayleighScattering;}
    G4OpMieHG* GetMieHGScatteringProcess() {return theMieHGScatteringProcess;}
    G4OpBoundaryProcess* GetBoundaryProcess() { return theBoundaryProcess;}

    void SetNbOfPhotonsCerenkov(G4int);

private:

    G4OpWLS*             theWLSProcess;
    G4Cerenkov*          theCerenkovProcess;
    G4Scintillation*     theScintProcess;
    G4OpAbsorption*      theAbsorptionProcess;
    G4OpRayleigh*        theRayleighScattering;
    G4OpMieHG*           theMieHGScatteringProcess;
    G4OpBoundaryProcess* theBoundaryProcess;
 
    G4bool AbsorptionOn;

};
#endif
