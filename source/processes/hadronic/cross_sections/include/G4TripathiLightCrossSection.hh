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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the	      *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 17191/03/NL/LvH (Aurora Programme). 		      *
// *                                                                  *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4TripathiLightCrossSection_h
#define G4TripathiLightCrossSection_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TripathiLightCrossSection.hh
//
// Version:		B.1
// Date:		15/04/04
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		17191/03/NL/LvH
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 6 October 2003, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// 15 March 2004, P R Truscott, QinetiQ Ltd, UK
// Beta release
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
// Implementation of formulas of Tripathi, Cucinotta and Wilson, NASA Technical
// Paper TP-1999-209726 to calculate cross-sections for nuclear-nuclear
// inelastic scattering for light nuclear systems.
//
// Class Description - End
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4WilsonRadius.hh"
///////////////////////////////////////////////////////////////////////////////
//
class G4TripathiLightCrossSection : public G4VCrossSectionDataSet
{
  public:
    G4TripathiLightCrossSection();
    ~G4TripathiLightCrossSection();
    virtual G4bool IsApplicable(const G4DynamicParticle* theProjectile,
      const G4Element* theTarget);

    virtual G4double GetCrossSection(const G4DynamicParticle* theProjectile,
      const G4Element* theTarget, G4double theTemperature);

    virtual void BuildPhysicsTable(const G4ParticleDefinition&)
    {}

    virtual void DumpPhysicsTable(const G4ParticleDefinition&)
    {G4cout << "G4TripathiLightCrossSection: uses formula"<<G4endl;}
    void SetLowEnergyCheck (G4bool);

  private:
    G4WilsonRadius *theWilsonRadius;
    G4double       r_0;
    G4double       third;
    G4bool         lowEnergyCheck;
};
///////////////////////////////////////////////////////////////////////////////
//
#endif
