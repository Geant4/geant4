// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Ions.hh,v 1.6 1999-12-15 14:51:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 27 June 1998
// ----------------------------------------------------------------
//      Add excitation energy         17 Aug. 1999 H.Kurashige
//


#ifndef G4Ions_h
#define G4Ions_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleWithCuts.hh"

// ######################################################################
// ###                          Ions                                 ###
// ######################################################################

class G4Ions : public G4ParticleWithCuts
{
 // Class Description
 //  All nuclei/ions created on the fly are objects of this class
 //  This class has Excitation Energy in addition to the normal particle
 //

 private:
   G4double  theIonsLengthCut;
   G4double* theIonsKineticEnergyCuts;

 public: //With Description
   G4Ions(
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
   virtual    			~G4Ions(){};
   G4Ions*    			IonsDefinition();
   G4Ions*    			Ions();

 public: //With Description
   virtual G4double 	   	GetCuts() {return theIonsLengthCut;}   
   virtual const G4double* 	GetCutsInEnergy() {return theIonsKineticEnergyCuts;};

   virtual void 		SetCuts(G4double aCut); 


 public:  //With Description
   G4int    GetAtomicNumber() const;
   G4int    GetAtomicMass() const;

   G4double GetExcitationEnergy() const ; 
   void     SetExcitationEnergy(G4double value);
  
  private:
   G4double theExcitationEnergy; 

};

inline
 G4Ions* G4Ions::Ions() 
{
  return this;
}

inline
 G4int G4Ions::GetAtomicNumber() const 
{
  return int(GetPDGCharge()/eplus); 
}

inline
 G4int G4Ions::GetAtomicMass() const 
{
  return GetBaryonNumber();
}

inline
 G4double G4Ions::GetExcitationEnergy() const 
{
  return theExcitationEnergy;
}

inline
 void G4Ions::SetExcitationEnergy(G4double value) 
{
  theExcitationEnergy = value;
}

inline 
 void G4Ions::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theIonsLengthCut = theCutInMaxInteractionLength;  
  theIonsKineticEnergyCuts = theKineticEnergyCuts;
}

#endif






