// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPElasticFS.hh,v 1.1 1999-01-07 16:12:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPElasticFS_h
#define G4NeutronHPElasticFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPLegendreStore.hh"
#include "G4NeutronHPPartial.hh"
#include "G4NeutronHPFastLegendre.hh"
#include "G4NeutronHPInterpolator.hh"
#include "G4NeutronHPNames.hh"

class G4NeutronHPElasticFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPElasticFS()
  {
    hasXsec = false; 
    theCoefficients = NULL;
    theProbArray = NULL;
  }
  ~G4NeutronHPElasticFS()
  {
    if(theCoefficients!=NULL) delete theCoefficients;
    if(theProbArray!=NULL) delete theProbArray;
  }
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4ParticleChange * ApplyYourself(const G4Track & theTrack);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPElasticFS * theNew = new G4NeutronHPElasticFS;
   return theNew;
  }
  
  private:
  G4int repFlag;    // Legendre coeff(1), or probability array(2), or isotropic(0).
  G4double targetMass; // in neutronmass units.
  G4int frameFlag;  // CMS or Lab system.
  
  G4NeutronHPLegendreStore * theCoefficients; // the legendre coefficients
  G4NeutronHPPartial * theProbArray; // the probability array p,costh for energy
  G4NeutronHPInterpolator theInt; // interpolation
  
  G4NeutronHPFastLegendre theLegend; // fast look-up for leg-integrals
};
#endif
