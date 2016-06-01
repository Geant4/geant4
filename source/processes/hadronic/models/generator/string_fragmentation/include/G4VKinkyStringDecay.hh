// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VKinkyStringDecay.hh,v 1.2 1998/11/03 10:10:34 maxim Exp $
// GEANT4 tag $Name: geant4-00 $
//  Maxim Komogorov
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, Maxim Komogorov, 10-Oct-1998
// -----------------------------------------------------------------------------

#ifndef G4VKinkyStringDecay_h
#define G4VKinkyStringDecay_h 1

#include "G4VLongitudinalStringDecay.hh"

//*********************************************************************************************** 

class G4VKinkyStringDecay 
    {

// Constructors   
public:
    G4VKinkyStringDecay(G4VLongitudinalStringDecay* theModal);
   ~G4VKinkyStringDecay() {};

// 
public:
    G4KineticTrackVector* FragmentString(const G4ExcitedString& String);
    virtual G4double GetLightConeGluonZ(G4double zmin, G4double zmax);
    void SetLongitudinalStringDecay(G4VLongitudinalStringDecay*);

private:
   G4VLongitudinalStringDecay* theLongitudinalStringDecay;  
   
   };

//*****************************************************************************************

inline void G4VKinkyStringDecay::SetLongitudinalStringDecay(G4VLongitudinalStringDecay* theModal)
   {
   theLongitudinalStringDecay = theModal;
   }

//*****************************************************************************************************

#endif
