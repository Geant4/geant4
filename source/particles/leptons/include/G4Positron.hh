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
// $Id: G4Positron.hh,v 1.7 2001-10-16 08:16:16 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      4-th April 1996, G.Cosmo
// ****************************************************************
//  Added particle definitions, H.Kurashige, 19 April 1996
//  Revised, G.Cosmo, 6 June 1996
//  Added SetCuts, L.Urban, 12 June 1996
//  Revised, Hisaya Kurashige, 4 July 1996
//  Revised, Hisaya Kurashige, 7 July 1996
//  Added not static GetEnergyCuts() and GetLengthCuts(), G.Cosmo, 11 July 1996
//  Revised, Hisaya Kurashige, 15 Dec 1996
// ----------------------------------------------------------------
// Each class inheriting from G4VLepton
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4Positron_h
#define G4Positron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VLepton.hh"
#include "G4Electron.hh"

// ######################################################################
// ###                         POSITRON                               ###
// ######################################################################

class G4Positron : public G4VLepton
{
 private:
   static G4Positron thePositron;
 
 protected:  
   G4double ComputeLoss(G4double AtomicNumber, G4double KineticEnergy) const;

   void BuildRangeVector( const G4Material* aMaterial,
			  const G4LossTable* aLossTable,
			  G4double       maxEnergy,     
			  G4double       aMass,
                          G4RangeVector* rangeVector
                         );

   friend void G4Electron::BuildRangeVector(
				  const G4Material* aMaterial,
				  const G4LossTable* aLossTable,
				  G4double       maxEnergy,     
				  G4double       aMass,
                                  G4RangeVector* rangeVector
                         );


 private:   //hide constructor as private 
   G4Positron(	
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable
   );

 public:
   virtual ~G4Positron(){};	
 
   static G4Positron* PositronDefinition();
   static G4Positron* Positron();
};


#endif



