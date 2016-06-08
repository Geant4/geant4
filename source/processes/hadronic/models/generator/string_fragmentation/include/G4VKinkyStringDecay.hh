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
//
// $Id: G4VKinkyStringDecay.hh,v 1.5 2001/10/05 16:17:36 hpw Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//  Maxim Komogorov
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
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
   virtual ~G4VKinkyStringDecay() {};

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
