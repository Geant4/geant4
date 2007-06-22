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
// $Id: G4NeutronHPAngular.hh,v 1.12 2007-06-22 09:23:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 070613 fix memory leaking by T. Koi
//
#ifndef G4NeutronHPAngular_h
#define G4NeutronHPAngular_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
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
// TKDB
      theCoefficients = 0;
      theProbArray = 0;
  } 
  ~G4NeutronHPAngular()
   {
// TKDB
      delete theCoefficients;
      delete theProbArray;
   }
  
  void Init(std::ifstream & aDataFile);
  
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
