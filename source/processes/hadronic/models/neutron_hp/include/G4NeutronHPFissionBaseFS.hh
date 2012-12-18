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
//
#ifndef G4NeutronHPFissionBaseFS_h
#define G4NeutronHPFissionBaseFS_h 1

#include "globals.hh"
#include "G4ReactionProduct.hh"
#include "G4DynamicParticleVector.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPVector.hh"
#include "G4NeutronHPEnergyDistribution.hh"
#include "G4NeutronHPAngular.hh"

class G4NeutronHPFissionBaseFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPFissionBaseFS()
  { 
    hasXsec = true; 
    theXsection = new G4NeutronHPVector;
  }
  virtual ~G4NeutronHPFissionBaseFS()
  {
    delete theXsection;
  }

  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & bit);

  G4DynamicParticleVector * ApplyYourself(G4int Prompt);

  virtual G4double GetXsec(G4double anEnergy)
  {
    return std::max(0., theXsection->GetY(anEnergy));
  }
  virtual G4NeutronHPVector * GetXsec() { return theXsection; }

  inline void SetNeutron(const G4ReactionProduct & aNeutron)
                        { 
                          theNeutron = aNeutron;
                          theAngularDistribution.SetNeutron(aNeutron);
                        }
  
  inline void SetTarget(const G4ReactionProduct & aTarget)
                        { 
                          theTarget = aTarget; 
                          theAngularDistribution.SetTarget(aTarget);
                        }
  
  private:
  
  G4HadFinalState * ApplyYourself(const G4HadProjectile & ) {return 0;}
  
  G4NeutronHPVector * theXsection;
  G4NeutronHPEnergyDistribution theEnergyDistribution;
  G4NeutronHPAngular theAngularDistribution;
  
  G4ReactionProduct theNeutron;
  G4ReactionProduct theTarget;

  private:
  
};
#endif
