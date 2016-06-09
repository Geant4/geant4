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
// $Id: G4SDParticleWithEnergyFilter.hh,v 1.2 2005/11/17 22:53:38 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#ifndef G4SDParticleWithEnergyFilter_h
#define G4SDParticleWithEnergyFilter_h 1

class G4Step;
class G4ParticleDefinition;
#include "globals.hh"
#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDKineticEnergyFilter.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector. 
//  This class filters steps by partilce definition and kinetic energy.
//
// Created: 2005-11-14  Tsukasa ASO.
// 
///////////////////////////////////////////////////////////////////////////////

class G4SDParticleWithEnergyFilter : public G4VSDFilter 
{
  public: // with description
      G4SDParticleWithEnergyFilter(G4String name,
				   G4double elow=0.0, G4double ehigh=DBL_MAX);
      virtual ~G4SDParticleWithEnergyFilter();

  public: // with description
      virtual G4bool Accept(const G4Step*) const;

      void add(const G4String& particleName);
      // add the particle into accepatable particle list.
      //
      void SetKineticEnergy(G4double elow, G4double ehigh);
      // Set acceptable kinetic energy range.
      //
      void show();

  private:
     G4SDParticleFilter* fParticleFilter;
     G4SDKineticEnergyFilter* fKineticFilter;
};

#endif

