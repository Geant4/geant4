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
// Short description: Interface of QGSC to CHIPS (all soft hadrons) 
// ----------------------------------------------------------------
//
#ifndef G4StringChipsParticleLevelInterface_h
#define G4StringChipsParticleLevelInterface_h

#include "G4VIntraNuclearTransportModel.hh"
#include "G4ChiralInvariantPhaseSpace.hh"

// Open this if you wish histogramming of the impact parameter issue
//#define hdebug_SCPLI

class G4StringChipsParticleLevelInterface : public G4VIntraNuclearTransportModel
{
  public:
    G4StringChipsParticleLevelInterface();
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

#ifdef hdebug_SCPLI
    //Static variables for histogramming
    static const G4int nbh;
    static       G4double bhmax;
    static       G4double bhdb;
    static       G4double ehmax;
    static       G4double ehde;
    static       G4double toth;
    static       G4int    bover;
    static       G4int    eover;
    static       G4int*   bhis;
    static       G4int*   ehis;
  public:
    //Static functions
    static void Reset()
    {
      bhdb=bhmax/nbh;
      ehde=bhmax/nbh;
      toth=0.;
      bover=0;
      eover=0;
      for(G4int i=0; i<nbh; i++)
      {
        bhis[i]=0;
        ehis[i]=0;
      }
    }
    static void SetMaxB(G4double mB) {bhmax=mB;}
    static void SetMaxE(G4double mE) {ehmax=mE;}
    static G4int GetB(G4int i){return bhis[i];}
    static G4int GetE(G4int i){return ehis[i];}
    static G4double GetTot()  {return toth;}
    static G4int GetNbn()     {return nbh;}
    static G4double GetDB()   {return bhdb;}
    static G4double GetDE()   {return ehde;}
    static G4int GetBov()     {return bover;}
    static G4int GetEov()     {return eover;}
#endif
};
#endif
