// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPAngular.hh,v 1.3 1999-07-02 09:58:23 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPAngular_h
#define G4NeutronHPAngular_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "G4ReactionProduct.hh"
#include "Randomize.hh"
#include "G4NeutronHPLegendreStore.hh"
#include "G4NeutronHPPartial.hh"

class G4NeutronHPAngular
{
    public:
    
  G4NeutronHPAngular()
  {
    theAngularDistributionType = 0;
    theIsoFlag = false;
  } 
  ~G4NeutronHPAngular(){}
  
  void Init(ifstream & aDataFile);
  
  void SampleAndUpdate(G4ReactionProduct & aNeutron);
    
  void SetTarget(const G4ReactionProduct & aTarget) { theTarget = aTarget; }

  void SetNeutron(const G4ReactionProduct & aNeutron) { theNeutron = aNeutron; }

  inline G4double GetTargetMass() { return targetMass; }

  private:
  
  // the type of distribution; currently 
  // isotropic (0), 
  // and legendre representation (1)
  // probability distribution (2)
  // are supported
  
  G4int theAngularDistributionType;
  G4int frameFlag; // 1=Lab, 2=CMS
    
  G4bool theIsoFlag; // isotropic or not?
  
  G4NeutronHPLegendreStore * theCoefficients; // the legendre coefficients

  G4NeutronHPPartial * theProbArray; // the probability array p,costh for energy

  private:
  
  G4double targetMass;

  G4ReactionProduct theTarget;
  G4ReactionProduct theNeutron;
};

#endif
