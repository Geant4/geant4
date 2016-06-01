// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VKinkyStringDecay.cc,v 1.4 1998/12/13 15:46:21 pia Exp $
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

#include "G4VKinkyStringDecay.hh"

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
    while (G4UniformRand() > yf); 
    return z;
    }

//*****************************************************************************************************

G4KineticTrackVector* G4VKinkyStringDecay::FragmentString(const G4ExcitedString& String) 
    {
    G4LorentzVector Mom = String.GetGluon()->Get4Momentum();
    G4ThreeVector Pos = String.GetGluon()->GetPosition();
    G4int QuarkEncoding = theLongitudinalStringDecay->SampleQuarkFlavor();
    G4double Px, Py; 
    theLongitudinalStringDecay->SampleQuarkPt(&Px, &Py);
    G4double Pt2 = Px*Px + Py*Py;
    G4double z = GetLightConeGluonZ(0, 1);
    G4double w = Mom.e() + Mom.pz();
    //... now compute quark longitudinal momentum and energy

    G4double Pz = (z*w - Pt2/(z*w))*0.5;
    G4double E  = (z*w + Pt2/(z*w))*0.5;
    
    G4Parton* AntiColor = new G4Parton(-QuarkEncoding);
    AntiColor->SetPosition(Pos);
    G4LorentzVector AntiColorMom(-Px, -Py, -Pz, E);
    AntiColor->Set4Momentum(AntiColorMom);
    G4Parton* Color = new G4Parton(*String.GetColorParton());
    G4ExcitedString Str1(Color, AntiColor, String.GetDirection());
    G4KineticTrackVector* KTV1 = theLongitudinalStringDecay->FragmentString(Str1);

    Color = new G4Parton(QuarkEncoding);
    Color->SetPosition(Pos);
    G4LorentzVector ColorMom(Px, Py, Pz, E);
    Color->Set4Momentum(ColorMom);
    AntiColor = new G4Parton(*String.GetAntiColorParton());
    G4ExcitedString Str2(Color, AntiColor, String.GetDirection());
    G4KineticTrackVector* KTV2 = theLongitudinalStringDecay->FragmentString(Str2);

    if (KTV1 && KTV2)
         while(!KTV2->isEmpty())
             KTV1->insert(KTV2->removeLast());
    return KTV1;            
    } 
 
//*****************************************************************************************************
 
 
 
 
 
 
 
 
 
 
 
 
 
 
