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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NeutronHPInelasticBaseFS.hh,v 1.6 2001-07-26 09:28:05 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPInelasticBaseFS_h
#define G4NeutronHPInelasticBaseFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPAngular.hh"
#include "G4NeutronHPEnergyDistribution.hh"
#include "G4NeutronHPEnAngCorrelation.hh"
#include "G4NeutronHPPhotonDist.hh"
#include "G4NeutronHPDeExGammas.hh"

class G4NeutronHPInelasticBaseFS : public G4NeutronHPFinalState
{
  public:
    
  G4NeutronHPInelasticBaseFS()
  {
    hasXsec = true; 
    theXsection = new G4NeutronHPVector;
    
    theEnergyDistribution = NULL;
    theFinalStatePhotons = NULL;
    theEnergyAngData = NULL;
    theAngularDistribution = NULL;
  }
  virtual ~G4NeutronHPInelasticBaseFS()
  {
    delete theXsection;
    if(theEnergyDistribution!=NULL) delete theEnergyDistribution;
    if(theFinalStatePhotons!=NULL) delete theFinalStatePhotons;
    if(theEnergyAngData!=NULL) delete theEnergyAngData;
    if(theAngularDistribution!=NULL) delete theAngularDistribution;
  }
  
  void Init (G4double A, G4double Z, G4String & dirName, G4String & bit);
  void BaseApply(const G4Track & theTrack, G4ParticleDefinition ** theDefs, G4int nDef);
  void InitGammas(G4double AR, G4double ZR);
  virtual G4ParticleChange * ApplyYourself(const G4Track & theTrack) = 0;
  virtual G4NeutronHPFinalState * New() = 0;
  
  virtual G4double GetXsec(G4double anEnergy)
  {
    return G4std::max(0., theXsection->GetY(anEnergy));
  }
  virtual G4NeutronHPVector * GetXsec() { return theXsection; }

  protected:
  
  G4NeutronHPVector * theXsection;
  G4NeutronHPEnergyDistribution * theEnergyDistribution;
  G4NeutronHPAngular * theAngularDistribution;
  G4NeutronHPEnAngCorrelation * theEnergyAngData;
  
  G4NeutronHPPhotonDist * theFinalStatePhotons;
  G4double theNuclearMassDifference;
  G4NeutronHPDeExGammas theGammas;
  G4String gammaPath;

};
#endif
