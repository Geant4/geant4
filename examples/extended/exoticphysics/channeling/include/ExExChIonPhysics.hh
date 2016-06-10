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


#ifndef ExExChIonPhysics_h
#define ExExChIonPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4HadronicInteraction;
class G4VCrossSectionDataSet;
class G4VComponentCrossSection;
class G4FTFBuilder;
class G4BinaryLightIonReaction;

class ExExChIonPhysics : public G4VPhysicsConstructor
{
public:
    
    ExExChIonPhysics(G4int ver = 0);
    ExExChIonPhysics(const G4String& nname);
    virtual ~ExExChIonPhysics();
    
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    void ConstructParticle();
    void ConstructProcess();
    
private:
    
    void AddProcess(const G4String&, G4ParticleDefinition*, G4bool isIon);
    
    static G4ThreadLocal G4VCrossSectionDataSet* theNuclNuclData;
    static G4ThreadLocal G4VComponentCrossSection* theGGNuclNuclXS;
    
    static G4ThreadLocal G4BinaryLightIonReaction* theIonBC;
    static G4ThreadLocal G4HadronicInteraction*    theFTFP;
    static G4ThreadLocal G4FTFBuilder*             theBuilder;
    
    G4int  verbose;
    static G4ThreadLocal G4bool wasActivated;
};


#endif








