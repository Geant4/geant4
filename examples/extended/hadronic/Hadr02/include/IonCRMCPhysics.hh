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
/// \file hadronic/Hadr02/include/IonCRMCPhysics.hh
/// \brief Definition of the IonCRMCPhysics class
//
//
//---------------------------------------------------------------------------
//
// Header:    IonCRMCPhysics
//
// Author:    2018 Alberto Ribon
//
// This is a variant of G4IonPhysics whereby CRMC is used for modeling
// final-state for ion nuclear inelastic interactions at very high energies.
// The inelastic hadronic cross sections are the same as in G4IonPhysics.
//
// Modified:     
//
// ------------------------------------------------------------
//
#ifndef G4IonCRMCPhysics_h
#define G4IonCRMCPhysics_h 1

#include "G4VHadronPhysics.hh"
#include "globals.hh"

class G4HadronicInteraction;
class G4VCrossSectionDataSet;
class G4VComponentCrossSection;
class G4FTFBuilder;
class G4BinaryLightIonReaction;
class G4CRMCModel;


class IonCRMCPhysics : public G4VPhysicsConstructor {
  public:
    IonCRMCPhysics( G4int ver = 0 );
    virtual ~IonCRMCPhysics();
    void ConstructParticle();
    void ConstructProcess();
  private:
    void AddProcess( const G4String& , G4ParticleDefinition* , G4bool isIon );
    static G4ThreadLocal G4VCrossSectionDataSet*   theNuclNuclData; 
    static G4ThreadLocal G4VComponentCrossSection* theGGNuclNuclXS;
    static G4ThreadLocal G4BinaryLightIonReaction* theIonBC;
    static G4ThreadLocal G4HadronicInteraction*    theFTFP;
    static G4ThreadLocal G4FTFBuilder*             theBuilder;
    static G4ThreadLocal G4CRMCModel*              theCRMC; 
    G4int  verbose;
    static G4ThreadLocal G4bool wasActivated;
};

#endif

