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
// $Id: G4XicPlus.hh,v 1.5 2001-07-11 10:01:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      4-th April 1996, G.Cosmo
// ****************************************************************
//  Added particle definitions, H.Kurashige, 14 Feb 19
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ----------------------------------------------------------------

// Each class inheriting from G4VBaryon
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4XicPlus_h
#define G4XicPlus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBaryon.hh"

// ######################################################################
// ###                          XicPlus                               ###
// ######################################################################

class G4XicPlus : public G4VBaryon
{
 private:
   static G4XicPlus theXicPlus;
   static G4double  theXicPlusLengthCut;
   static G4double* theXicPlusKineticEnergyCuts;

 private:
   G4XicPlus(
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
   virtual ~G4XicPlus(){}

   static G4XicPlus* XicPlusDefinition();
   static G4XicPlus* XicPlus();
   static G4double GetCuts() {return theXicPlusLengthCut;}   
   static G4double* GetCutsInEnergy() {return theXicPlusKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};

inline void G4XicPlus::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theXicPlusLengthCut = theCutInMaxInteractionLength;  
  theXicPlusKineticEnergyCuts = theKineticEnergyCuts;
  
}


#endif
