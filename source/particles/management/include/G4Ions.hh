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
// $Id: G4Ions.hh,v 1.9 2001-10-15 09:58:30 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
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
 G4double G4Ions::GetExcitationEnergy() const 
{
  return theExcitationEnergy;
}

inline
 void G4Ions::SetExcitationEnergy(G4double value) 
{
  theExcitationEnergy = value;
}

#endif
