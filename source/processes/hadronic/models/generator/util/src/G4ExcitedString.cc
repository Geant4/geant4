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
// $Id: G4ExcitedString.cc,v 1.6 2001/10/04 20:00:33 hpw Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4ExcitedString ----------------
//             by Gunter Folger, June 1998.
//       class for an excited string used by Parton String Models
// ------------------------------------------------------------


// G4ExcitedString
#include "G4ExcitedString.hh"
#include "g4std/algorithm"

//G4ExcitedString::G4ExcitedString(const G4ExcitedString &right)
//{}

G4ExcitedString::G4ExcitedString(G4Parton* Color, G4Parton* AntiColor, G4int Direction)
    {
    thePartons.push_back(Color);
    thePartons.push_back(AntiColor);
    thePosition = Color->GetPosition();
    theDirection = Direction;
    theTrack=0;
    }

G4ExcitedString::G4ExcitedString(G4Parton* Color, G4Parton* Gluon,  G4Parton* AntiColor, G4int Direction)
    {
    thePartons.push_back(Color);
    thePartons.push_back(Gluon);
    thePartons.push_back(AntiColor);
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
  G4std::for_each(thePartons.begin(), thePartons.end(), DeleteParton());
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
    for(unsigned int cParton = 0; cParton < thePartons.size() ; cParton++ )
        {
        G4LorentzVector Mom = thePartons[cParton]->Get4Momentum();
        Mom.boost(Velocity);
        thePartons[cParton]->Set4Momentum(Mom);
        }
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetColorParton(void) const
    {
    G4Parton * start = *(thePartons.begin());
    G4Parton * end = *(thePartons.end()-1);
    G4int Encoding = start->GetPDGcode();
    if (Encoding < -1000 || ((Encoding  < 1000) && (Encoding > 0)))
        return start;
    return end; 
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetGluon(void) const
    {
    return thePartons[1]; 
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetGluon(G4int GluonPos) const
    {
    return thePartons[1 + GluonPos]; 
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetAntiColorParton(void) const
    {
    G4Parton * start = *(thePartons.begin());
    G4Parton * end = *(thePartons.end()-1);
    G4int Encoding = start->GetPDGcode();
    if (Encoding < -1000 || ((Encoding  < 1000) && (Encoding > 0)))
        return end; 
    return start; 
    }

//---------------------------------------------------------------------------------

G4bool G4ExcitedString::IsItKinkyString(void) const
    {
    return (thePartons.size() > 2);    
    }

//---------------------------------------------------------------------------------

G4int G4ExcitedString::GetDirection(void) const
    {
    return theDirection;    
    }

//*********************************************************************************

G4Parton* G4ExcitedString::GetLeftParton(void) const
    {
    return *thePartons.begin(); 
    }

//---------------------------------------------------------------------------------

G4Parton* G4ExcitedString::GetRightParton(void) const
    {
    return *(thePartons.end()-1); 
    }

//*********************************************************************************
