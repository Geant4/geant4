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
// $Id: G4NeutronHPInelasticCompFS.hh,v 1.6 2001-07-11 10:06:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPInelasticCompFS_h
#define G4NeutronHPInelasticCompFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPAngular.hh"
#include "G4NeutronHPEnergyDistribution.hh"
#include "G4NeutronHPEnAngCorrelation.hh"
#include "G4NeutronHPPhotonDist.hh"
#include "G4NeutronHPDeExGammas.hh"

class G4NeutronHPInelasticCompFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPInelasticCompFS()
  {
    for(G4int i=0; i<51; i++)
    {
      hasXsec = true; 
      theXsection[i] = NULL;
      theEnergyDistribution[i] = NULL;
      theAngularDistribution[i] = NULL;
      theEnergyAngData[i] = NULL;
      theFinalStatePhotons[i] = NULL;
    }
  }
  virtual ~G4NeutronHPInelasticCompFS()
  {
    for(G4int i=0; i<51; i++)
    {
      if(theXsection[i] != NULL) delete theXsection[i];
      if(theEnergyDistribution[i] != NULL) delete theEnergyDistribution[i];
      if(theAngularDistribution[i] != NULL) delete theAngularDistribution[i];
      if(theEnergyAngData[i] != NULL) delete theEnergyAngData[i];
      if(theFinalStatePhotons[i] != NULL) delete theFinalStatePhotons[i];
    }
  }
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aSFType);
  void InitGammas(G4double AR, G4double ZR);
  virtual G4ParticleChange * ApplyYourself(const G4Track & theTrack) = 0;
  virtual G4NeutronHPFinalState * New() = 0;
  virtual G4double GetXsec(G4double anEnergy)
  {
    return G4std::max(0., theXsection[50]->GetY(anEnergy));
  }
  virtual G4NeutronHPVector * GetXsec() { return theXsection[50]; }
  G4int SelectExitChannel(G4double eKinetic);
  void CompositeApply(const G4Track & theTrack, G4ParticleDefinition * aHadron);
  inline void InitDistributionInitialState(G4ReactionProduct & aNeutron, 
                                           G4ReactionProduct & aTarget, 
                                           G4int it)
  {
    if(theAngularDistribution[it]!=NULL) 
    {
      theAngularDistribution[it]->SetTarget(aTarget);
      theAngularDistribution[it]->SetNeutron(aNeutron);
    }
    if(theEnergyAngData[it]!=NULL)
    {
      theEnergyAngData[it]->SetTarget(aTarget);
      theEnergyAngData[it]->SetNeutron(aNeutron);
    }
  }
  
  protected:
  
  G4NeutronHPVector * theXsection[51];
  G4NeutronHPEnergyDistribution * theEnergyDistribution[51];
  G4NeutronHPAngular * theAngularDistribution[51];
  G4NeutronHPEnAngCorrelation * theEnergyAngData[51];
  
  G4NeutronHPPhotonDist * theFinalStatePhotons[51];
  
  G4NeutronHPDeExGammas theGammas;
  G4String gammaPath;
  
  G4double theCurrentA;
  G4double theCurrentZ;
};
#endif
