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
// $Id: G4GenericIon.hh,v 1.7 2001-10-16 08:16:05 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      4-th Dec 1998, H.Kurashige
// ****************************************************************
// This class is used only by G4IonTable and not for tracking
// G4IonTable creates various ions other than alpha,deuteron,triton,
// and He3. Processes for these ions will be same as ones for 
// this "GenericIon". So, user should register processes for ions
// to this class in his/her UserPhysicsList

#ifndef G4GenericIon_h
#define G4GenericIon_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIon.hh"

// ######################################################################
// ###                          GenericIon                            ###
// ######################################################################

class G4GenericIon : public G4VIon
{
 private:
   static G4GenericIon theGenericIon;

 public:
   G4GenericIon(
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
   virtual ~G4GenericIon();

   static G4GenericIon*    GenericIonDefinition();
   static G4GenericIon*    GenericIon();
};

#endif
