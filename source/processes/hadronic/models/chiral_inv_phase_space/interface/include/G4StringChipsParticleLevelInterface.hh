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
#ifndef G4StringChipsParticleLevelInterface_h
#define G4StringChipsParticleLevelInterface_h

#include "G4VIntraNuclearTransportModel.hh"
#include "G4ChiralInvariantPhaseSpace.hh"

class G4StringChipsParticleLevelInterface : public G4VIntraNuclearTransportModel
{
  public:
    G4StringChipsParticleLevelInterface();
    virtual G4VParticleChange* ApplyYourself(const G4Track& aTrack, 
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
