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
// $Id: G4He3.hh,v 1.5 2001-07-11 10:01:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      24-th April 1998,H.Kurashige
// ****************************************************************
// ----------------------------------------------------------------------
// Each class inheriting from G4VIon
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4He3_h
#define G4He3_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIon.hh"

// ######################################################################
// ###                          He3                                 ###
// ######################################################################

class G4He3 : public G4VIon
{
 private:
   static G4He3 theHe3;
   static G4double  theHe3LengthCut;
   static G4double* theHe3KineticEnergyCuts;

 public:
   G4He3(
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
   virtual ~G4He3();

   static G4He3*    He3Definition();
   static G4He3*    He3() {return &theHe3;}
   static G4double GetCuts() {return theHe3LengthCut;}   
   static G4double* GetCutsInEnergy() {return theHe3KineticEnergyCuts;};

   void SetCuts(G4double aCut); 
   virtual void RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy );
};

inline void G4He3::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theHe3LengthCut = theCutInMaxInteractionLength;  
  theHe3KineticEnergyCuts = theKineticEnergyCuts;
  
}
inline void G4He3::RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy )
{
  G4ParticleWithCuts::RestoreCuts(cutInLength, cutInEnergy);
  theHe3LengthCut = theCutInMaxInteractionLength;  
  theHe3KineticEnergyCuts = theKineticEnergyCuts;
}
#endif

