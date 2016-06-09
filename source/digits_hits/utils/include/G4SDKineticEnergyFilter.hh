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
// $Id: G4SDKineticEnergyFilter.hh,v 1.2 2005/11/17 22:53:38 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#ifndef G4SDKineticEnergyFilter_h
#define G4SDKineticEnergyFilter_h 1

#include "globals.hh"
#include "G4VSDFilter.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector. 
//
//  This filter accepts particles defined energy range.
//  The energy range is given at constructor, or Set methods.
//
//
//
// Created: 2005-11-14  Tsukasa ASO.
// 
///////////////////////////////////////////////////////////////////////////////

class G4SDKineticEnergyFilter : public G4VSDFilter 
{

//-------
  public: // with description
      G4SDKineticEnergyFilter(G4String name,
			      G4double elow=0.0, 
			      G4double ehigh=DBL_MAX);
      // Constructor. Filter name and kinetic energy range( elow, ehigh).

     virtual ~G4SDKineticEnergyFilter();

  public: // with description
     virtual G4bool Accept(const G4Step*) const;

     void SetKineticEnergy(G4double elow, G4double ehigh);
     void SetLowEnergy(G4double elow);
     void SetHighEnergy(G4double ehigh);
     // Set methods for kinetic energy range.
     //
     void show();

  private:
     G4double fLowEnergy;
     G4double fHighEnergy;

};

#endif

