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
// $Id: G4DsMesonMinus.hh,v 1.4 2001-07-11 10:01:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
//      Created,             Hisaya Kurashige, 15 June 1997
// **********************************************************************
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ----------------------------------------------------------------

// Each class inheriting from G4VMeson
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4DsMesonMinus_h
#define G4DsMesonMinus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VMeson.hh"

// ######################################################################
// ###                       DsMesonMinus                             ###
// ######################################################################

class G4DsMesonMinus : public G4VMeson
{
 private:
   static G4DsMesonMinus theDsMesonMinus;
   static G4double  theDsMesonMinusLengthCut;
   static G4double* theDsMesonMinusKineticEnergyCuts;

 private: // constructors are hide as private  
   G4DsMesonMinus(
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
  virtual ~G4DsMesonMinus(){}

   static G4DsMesonMinus* DsMesonMinusDefinition();
   static G4DsMesonMinus* DsMesonMinus();
   static G4double GetCuts() {return theDsMesonMinusLengthCut;}   
   static G4double* GetCutsInEnergy() {return theDsMesonMinusKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};

inline void G4DsMesonMinus::SetCuts(G4double aCut)
{
  G4ParticleWithCuts::SetCuts(aCut);
  theDsMesonMinusLengthCut = theCutInMaxInteractionLength;  
  theDsMesonMinusKineticEnergyCuts = theKineticEnergyCuts;
}

#endif
