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
// $Id: G4NeutronHPInelasticBaseFS.hh,v 1.12 2007-06-06 12:45:13 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPInelasticBaseFS_h
#define G4NeutronHPInelasticBaseFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
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
    
    theEnergyDistribution = 0;
    theFinalStatePhotons = 0;
    theEnergyAngData = 0;
    theAngularDistribution = 0;
  }
  virtual ~G4NeutronHPInelasticBaseFS()
  {
    delete theXsection;
    if(theEnergyDistribution!=0) delete theEnergyDistribution;
    if(theFinalStatePhotons!=0) delete theFinalStatePhotons;
    if(theEnergyAngData!=0) delete theEnergyAngData;
    if(theAngularDistribution!=0) delete theAngularDistribution;
  }
  
  void Init (G4double A, G4double Z, G4String & dirName, G4String & bit);
  void BaseApply(const G4HadProjectile & theTrack, G4ParticleDefinition ** theDefs, G4int nDef);
  void InitGammas(G4double AR, G4double ZR);
  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack) = 0;
  virtual G4NeutronHPFinalState * New() = 0;
  
  virtual G4double GetXsec(G4double anEnergy)
  {
    return std::max(0., theXsection->GetY(anEnergy));
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
