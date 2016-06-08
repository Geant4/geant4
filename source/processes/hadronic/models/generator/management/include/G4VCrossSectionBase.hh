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
#ifndef G4VCrossSectionBase_h
#define G4VCrossSectionBase_h 1

#include "globals.hh"

const G4int maxpsig = 13;
const G4int maxreac = 14;

//****************************************************************************************

class G4VCrossSectionBase 
   {
public:
   //Constructors
   G4VCrossSectionBase();
  virtual ~G4VCrossSectionBase();

public:   
   virtual G4double GetTotCrossSection(G4int ProjectileEncoding, G4int TargetEncoding, G4double Energy)=0;
   virtual G4double GetElCrossSection (G4int ProjectileEncoding, G4int TargetEncoding, G4double Energy)=0;
   virtual G4double GetInCrossSection (G4int ProjectileEncoding, G4int TargetEncoding, G4double Energy)=0;   
   virtual G4double GetAnnihCrossSection(G4int ProjectileEncoding, G4int TargetEncoding, G4double
   Energy) = 0;
   };

//****************************************************************************************
#endif



