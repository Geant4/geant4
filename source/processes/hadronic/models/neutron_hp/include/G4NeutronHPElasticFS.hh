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
// $Id: G4NeutronHPElasticFS.hh,v 1.5 2001-07-26 09:27:57 hpw Exp $
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
