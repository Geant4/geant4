// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFissionBaseFS.hh,v 1.3 1999-07-02 09:59:02 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  void Init (G4double A, G4double Z, G4String & dirName, G4String & bit);

  G4DynamicParticleVector * ApplyYourself(G4int Prompt);

  virtual G4double GetXsec(G4double anEnergy)
  {
    return theXsection->GetY(anEnergy);
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
  
  G4ParticleChange * ApplyYourself(const G4Track & aTrack) {return NULL;}
  
  G4NeutronHPVector * theXsection;
  G4NeutronHPEnergyDistribution theEnergyDistribution;
  G4NeutronHPAngular theAngularDistribution;
  
  G4ReactionProduct theNeutron;
  G4ReactionProduct theTarget;

  private:
  
};
#endif
