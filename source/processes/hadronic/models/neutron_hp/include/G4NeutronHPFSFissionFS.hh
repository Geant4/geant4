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
// $Id: G4NeutronHPFSFissionFS.hh,v 1.11 2007-06-06 12:45:13 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFSFissionFS_h
#define G4NeutronHPFSFissionFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4DynamicParticleVector.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPNeutronYield.hh"
#include "G4NeutronHPVector.hh"
#include "G4NeutronHPFissionERelease.hh"
#include "G4NeutronHPEnergyDistribution.hh"
#include "G4NeutronHPPhotonDist.hh"
#include "G4NeutronHPAngular.hh"

class G4NeutronHPFSFissionFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPFSFissionFS(){ hasXsec = true; }
  ~G4NeutronHPFSFissionFS(){}
  
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  
  G4DynamicParticleVector * ApplyYourself(G4int Prompt, G4int delayed, G4double *decayconst);
  
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPFSFissionFS * theNew = new G4NeutronHPFSFissionFS;
   return theNew;
  }
  
  inline G4double GetMass(){ return targetMass; }
  
  void SampleNeutronMult(G4int&all, 
	  		 G4int&Prompt, 
			 G4int&delayed, 
			 G4double energy,
			 G4int off);
						 
  inline void SetNeutron(const G4ReactionProduct & aNeutron)
                        { 
                          theNeutron = aNeutron;
                          theNeutronAngularDis.SetNeutron(aNeutron);
                        }
  
  inline void SetTarget(const G4ReactionProduct & aTarget)
                        { 
                          theTarget = aTarget; 
                          theNeutronAngularDis.SetTarget(aTarget);
                        }
    
  G4DynamicParticleVector * GetPhotons();
  
  inline G4NeutronHPFissionERelease * GetEnergyRelease()
  {
    return &theEnergyRelease;
  }
  
  private:

  G4HadFinalState * ApplyYourself(const G4HadProjectile & ) { return 0; }  
  G4double targetMass;
  
  G4NeutronHPNeutronYield theFinalStateNeutrons;
  G4NeutronHPEnergyDistribution thePromptNeutronEnDis;
  G4NeutronHPEnergyDistribution theDelayedNeutronEnDis;
  G4NeutronHPAngular theNeutronAngularDis;
  
  G4NeutronHPPhotonDist theFinalStatePhotons;
  G4NeutronHPFissionERelease theEnergyRelease;
  
  G4ReactionProduct theNeutron;
  G4ReactionProduct theTarget;
  
  private:
  
  G4NeutronHPNames theNames;
  
};
#endif
