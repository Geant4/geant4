// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedString.cc,v 1.3 1999/12/15 14:52:51 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4ExcitedString ----------------
//             by Gunter Folger, June 1998.
//       class for an excited string used by Parton String Models
// ------------------------------------------------------------


// G4ExcitedString
#include "G4ExcitedString.hh"

//G4ExcitedString::G4ExcitedString(const G4ExcitedString &right)
//{}

G4ExcitedString::G4ExcitedString(G4Parton* Color, G4Parton* AntiColor, G4int Direction)
    {
    thePartons.insert(Color);
    thePartons.insert(AntiColor);
    thePosition = Color->GetPosition();
    theDirection = Direction;
    theTrack=0;
    }

G4ExcitedString::G4ExcitedString(G4Parton* Color, G4Parton* Gluon,  G4Parton* AntiColor, G4int Direction)
    {
    thePartons.insert(Color);
    thePartons.insert(Gluon);
    thePartons.insert(AntiColor);
    thePosition = Color->GetPosition();
    theDirection = Direction;
    theTrack=0;
    }

G4ExcitedString::G4ExcitedString(G4KineticTrack * track)
{
	thePosition = track->GetPosition();
	theTrack= track;
	theDirection=0;
}

G4ExcitedString::~G4ExcitedString()
{
	thePartons.clearAndDestroy();
}


//const G4ExcitedString & G4ExcitedString::operator=(const G4ExcitedString &right)
//{}


//int G4ExcitedString::operator==(const G4ExcitedString &right) const
//{}

//int G4ExcitedString::operator!=(const G4ExcitedString &right) const
//{}



// Additional Declarations


void G4ExcitedString::Boost(G4ThreeVector& Velocity)
    {
    for(G4int cParton = 0; cParton < thePartons.entries() ; cParton++ )
        {
        G4LorentzVector Mom = thePartons[cParton]->Get4Momentum();
        Mom.boost(Velocity);
        thePartons[cParton]->Set4Momentum(Mom);
        }
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetColorParton(void) const
    {
    G4int Encoding = thePartons.first()->GetPDGcode();
    if (Encoding < -1000 || ((Encoding  < 1000) && (Encoding > 0)))
        return thePartons.first(); 
    return thePartons.last(); 
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetGluon(void) const
    {
    return thePartons.at(1); 
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetGluon(G4int GluonPos) const
    {
    return thePartons.at(1 + GluonPos); 
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetAntiColorParton(void) const
    {
    G4int Encoding = thePartons.first()->GetPDGcode();
    if (Encoding < -1000 || ((Encoding  < 1000) && (Encoding > 0)))
        return thePartons.last(); 
    return thePartons.first(); 
    }

//---------------------------------------------------------------------------------

G4bool G4ExcitedString::IsItKinkyString(void) const
    {
    return (thePartons.entries() > 2);    
    }

//---------------------------------------------------------------------------------

G4int G4ExcitedString::GetDirection(void) const
    {
    return theDirection;    
    }

//*********************************************************************************

G4Parton* G4ExcitedString::GetLeftParton(void) const
    {
    return thePartons.first(); 
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetRightParton(void) const
    {
    return thePartons.last(); 
    }

//*********************************************************************************
