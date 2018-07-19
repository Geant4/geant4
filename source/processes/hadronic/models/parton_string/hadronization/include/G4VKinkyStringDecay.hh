//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VKinkyStringDecay.hh 102048 2016-12-19 09:02:38Z gcosmo $
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

//*****************************************************************************************

class G4VKinkyStringDecay 
{
  public:
    G4VKinkyStringDecay(G4VLongitudinalStringDecay* theModal);
    virtual ~G4VKinkyStringDecay() {};

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

//*****************************************************************************************

#endif

