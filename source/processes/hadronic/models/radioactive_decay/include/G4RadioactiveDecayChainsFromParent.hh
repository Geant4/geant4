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
#ifndef G4RadioactiveDecayChainsFromParent_h
#define G4RadioactiveDecayChainsFromParent_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              RadioactiveDecayRateVector.hh
//   renamed as         G4RadioactiveDecayChainsFromParent.hh  (D.H. wright  6 Oct 2017) 
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file     
//
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This class contains the decay times and coefficients for calculating      //
//  all the descendants in the decay chains of the named isotope.  These      //
//  data can be used to calculate their radioactivity at any given time.      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "globals.hh"
#include "G4RadioactiveDecayRatesToDaughter.hh"
#include <vector>

typedef std::vector<G4RadioactiveDecayRatesToDaughter> G4RadioactiveDecayRates;

class G4RadioactiveDecayChainsFromParent
{
  public:
    G4RadioactiveDecayChainsFromParent();
    virtual ~G4RadioactiveDecayChainsFromParent();
  
    G4RadioactiveDecayChainsFromParent(const G4RadioactiveDecayChainsFromParent&);
    G4RadioactiveDecayChainsFromParent& operator=(const G4RadioactiveDecayChainsFromParent&);
  
    // equality operators
    G4int operator==(const G4RadioactiveDecayChainsFromParent& right) const
      {return (this == &right);}
    G4int operator!=(const G4RadioactiveDecayChainsFromParent& right) const
      {return (this != &right);}
  
  public:

    inline G4String  GetIonName() const {return ionName;}
    inline void SetIonName(G4String name) {ionName = name;}

    // Retrieve the coefficients and decays of all descendants along the
    // decay chains
    inline G4RadioactiveDecayRates GetItsRates() const {return itsRates;}

    // Fill in the coefficients and decay times in the chains
    inline void SetItsRates(G4RadioactiveDecayRates arate) {itsRates = arate;}

protected:
    G4String ionName;
    G4RadioactiveDecayRates itsRates;

};
#endif

