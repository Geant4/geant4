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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPInelasticCompFS_h
#define G4ParticleHPInelasticCompFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPAngular.hh"
#include "G4ParticleHPEnergyDistribution.hh"
#include "G4ParticleHPEnAngCorrelation.hh"
#include "G4ParticleHPPhotonDist.hh"
#include "G4ParticleHPDeExGammas.hh"
#include "G4Nucleus.hh"

#include "G4NRESP71M03.hh"

class G4ParticleHPInelasticCompFS : public G4ParticleHPFinalState
{
  public:
  
  G4ParticleHPInelasticCompFS()
  {

    QI.resize(51);
    LR.resize(51);
    for(G4int i=0; i<51; i++)
    {
      hasXsec = true; 
      theXsection[i] = 0;
      theEnergyDistribution[i] = 0;
      theAngularDistribution[i] = 0;
      theEnergyAngData[i] = 0;
      theFinalStatePhotons[i] = 0;
      QI[i]=0.0;
      LR[i]=0;
    }

  }
  virtual ~G4ParticleHPInelasticCompFS()
  {
    for(G4int i=0; i<51; i++)
    {
      if(theXsection[i] != 0) delete theXsection[i];
      if(theEnergyDistribution[i] != 0) delete theEnergyDistribution[i];
      if(theAngularDistribution[i] != 0) delete theAngularDistribution[i];
      if(theEnergyAngData[i] != 0) delete theEnergyAngData[i];
      if(theFinalStatePhotons[i] != 0) delete theFinalStatePhotons[i];
    }
  }
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aSFType, G4ParticleDefinition*);
  void InitGammas(G4double AR, G4double ZR);
  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack) = 0;
  virtual G4ParticleHPFinalState * New() = 0;
  virtual G4double GetXsec(G4double anEnergy)
  {
    return std::max(0., theXsection[50]->GetY(anEnergy));
  }
  virtual G4ParticleHPVector * GetXsec() { return theXsection[50]; }
  G4int SelectExitChannel(G4double eKinetic);
  void CompositeApply(const G4HadProjectile & theTrack, G4ParticleDefinition * aHadron);
  inline void InitDistributionInitialState(G4ReactionProduct & anIncidentPart, 
                                           G4ReactionProduct & aTarget, 
                                           G4int it)
  {
    if(theAngularDistribution[it]!=0) 
    {
      theAngularDistribution[it]->SetTarget(aTarget);
      theAngularDistribution[it]->SetProjectileRP(anIncidentPart);
    }
    if(theEnergyAngData[it]!=0)
    {
      theEnergyAngData[it]->SetTarget(aTarget);
      theEnergyAngData[it]->SetProjectileRP(anIncidentPart);
    }
  }
  
protected:
  
  G4ParticleHPVector * theXsection[51];
  G4ParticleHPEnergyDistribution * theEnergyDistribution[51];
  G4ParticleHPAngular * theAngularDistribution[51];
  G4ParticleHPEnAngCorrelation * theEnergyAngData[51];
  
  G4ParticleHPPhotonDist * theFinalStatePhotons[51];
  
  G4ParticleHPDeExGammas theGammas;
  G4String gammaPath;
  
  //G4double theCurrentA;
  //G4double theCurrentZ;

   protected:
      std::vector < G4double >  QI;
      std::vector <G4int > LR;

   private:
      //                       proj                 targ                 had                  mu of had   
      void two_body_reaction ( G4DynamicParticle* , G4DynamicParticle* , G4DynamicParticle* , G4double mu ); 

      G4NRESP71M03 nresp71_model;
      G4bool use_nresp71_model( const G4ParticleDefinition* aDefinition, const G4int it , const G4ReactionProduct& theTarget , G4ReactionProduct& boosted);

};
#endif
