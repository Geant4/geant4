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
// $Id: G4VKinkyStringDecay.cc 102717 2017-02-20 10:37:13Z gcosmo $
//  Maxim Komogorov
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Oct-1998
// -----------------------------------------------------------------------------

#include "G4VKinkyStringDecay.hh"
#include "G4KineticTrackVector.hh"
#include "G4KineticTrack.hh"
#include "Randomize.hh"

//*****************************************************************************************************

G4VKinkyStringDecay::G4VKinkyStringDecay(G4VLongitudinalStringDecay* theModal)
   {
   this->SetLongitudinalStringDecay(theModal);
   }

//*****************************************************************************************************

G4double G4VKinkyStringDecay::GetLightConeGluonZ(G4double zmin, G4double zmax)
    {
    G4double z, yf;
    do {
       z = zmin + G4UniformRand()*(zmax-zmin);
       yf = z*z +sqr(1 - z);	
       } 
    while (G4UniformRand() > yf);  /* Loop checking, 07.08.2015, A.Ribon */ 
    return z;
    }

//*****************************************************************************************************

G4KineticTrackVector* G4VKinkyStringDecay::FragmentString(const G4ExcitedString& String) 
    {
    G4LorentzVector Mom = String.GetGluon()->Get4Momentum();
    G4ThreeVector Pos = String.GetGluon()->GetPosition();
    G4int QuarkEncoding = theLongitudinalStringDecay->SampleQuarkFlavor();
    G4ThreeVector Pquark=theLongitudinalStringDecay->SampleQuarkPt();
    G4double Pt2 = Pquark.mag2();
    G4double z = GetLightConeGluonZ(0, 1);
    G4double w = Mom.e() + Mom.pz();
    //... now compute quark longitudinal momentum and energy

    Pquark.setZ(  (z*w - Pt2/(z*w))*0.5);
    G4double E  = (z*w + Pt2/(z*w))*0.5;
    
    G4Parton* AntiColor = new G4Parton(-QuarkEncoding);
    AntiColor->SetPosition(Pos);
    G4LorentzVector AntiColorMom(-Pquark, E);
    AntiColor->Set4Momentum(AntiColorMom);
    G4Parton* Color = new G4Parton(*String.GetColorParton());
    G4ExcitedString Str1(Color, AntiColor, String.GetDirection());
    G4KineticTrackVector* KTV1 = theLongitudinalStringDecay->FragmentString(Str1);

    Color = new G4Parton(QuarkEncoding);
    Color->SetPosition(Pos);
    G4LorentzVector ColorMom(Pquark, E);
    Color->Set4Momentum(ColorMom);
    AntiColor = new G4Parton(*String.GetAntiColorParton());
    G4ExcitedString Str2(Color, AntiColor, String.GetDirection());
    G4KineticTrackVector* KTV2 = theLongitudinalStringDecay->FragmentString(Str2);

    if (KTV1 && KTV2)
         while(!KTV2->empty())  /* Loop checking, 07.08.2015, A.Ribon */
	 {
             KTV1->push_back(KTV2->back());
	     KTV1->erase(KTV1->end()-1);
	 }
    return KTV1;            
    } 
 
//*****************************************************************************************************
 
 
 
 
 
 
 
 
 
 
 
 
 
 
