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
// Short description: Interface of QGSC to CHIPS (Energy Flow of soft hadrons) 
// ---------------------------------------------------------------------------
//
#ifndef G4QStringChipsParticleLevelInterface_h
#define G4QStringChipsParticleLevelInterface_h

#include "G4VIntraNuclearTransportModel.hh"
#include "G4ChiralInvariantPhaseSpace.hh"

class G4QStringChipsParticleLevelInterface : public G4VIntraNuclearTransportModel
{
  public:
    G4QStringChipsParticleLevelInterface();
    virtual G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
                                             G4Nucleus& theNucleus);

    virtual G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
                                               G4V3DNucleus* theNucleus); 
  private:
  
    G4ChiralInvariantPhaseSpace theModel;
    G4double theEnergyLossPerFermi;
    
    G4double theInnerCoreDensityCut;
    
    G4double fractionOfSingleQuasiFreeNucleons;
    G4double fractionOfPairedQuasiFreeNucleons;
    G4double clusteringCoefficient;
    G4double temperature;
    G4double halfTheStrangenessOfSee;
    G4double etaToEtaPrime;
    G4double fusionToExchange;
    G4int nop;
};
#endif
