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
// $Id: G4NeutronHPAngular.hh,v 1.5 2001-07-11 10:06:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPAngular_h
#define G4NeutronHPAngular_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
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
  
  void Init(G4std::ifstream & aDataFile);
  
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
