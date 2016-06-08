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
// $Id: G4FragmentingString.cc,v 1.4 2002/06/13 09:04:12 jwellisc Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//


// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4FragmentingString ----------------
//             by Gunter Folger, September 2001.
//       class for an excited string used in Fragmention
// ------------------------------------------------------------


// G4FragmentingString
#include "G4FragmentingString.hh"
#include "G4ExcitedString.hh"

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------

G4FragmentingString::G4FragmentingString(const G4FragmentingString &old)
{
	LeftParton=old.LeftParton;
	RightParton=old.RightParton;
	Ptleft=old.Ptleft;
	Ptright=old.Ptright;
	Pplus=old.Pplus;
	Pminus=old.Pminus;
	decaying=old.decaying;
}

//---------------------------------------------------------------------------------

G4FragmentingString::G4FragmentingString(const G4ExcitedString &excited)
{
	LeftParton=excited.GetLeftParton()->GetDefinition();
	RightParton=excited.GetRightParton()->GetDefinition();
	Ptleft=excited.GetLeftParton()->Get4Momentum().vect();
	Ptleft.setZ(0.);
	Ptright=excited.GetRightParton()->Get4Momentum().vect();
	Ptright.setZ(0.);
	G4LorentzVector P=excited.Get4Momentum();
	Pplus =P.e() + P.pz();
	Pminus=P.e() - P.pz();
	decaying=None;
}

//---------------------------------------------------------------------------------

G4FragmentingString::G4FragmentingString(const G4FragmentingString &old,
					 G4ParticleDefinition * newdecay,
					 const G4LorentzVector *momentum)
{
	decaying=None;
	if ( old.decaying == Left )
	{
		RightParton= old.RightParton;
		Ptright    = old.Ptright;
		LeftParton = newdecay;
		Ptleft     = old.Ptleft - momentum->vect();
		Ptleft.setZ(0.);
	} else if ( old.decaying == Right )
	{
		RightParton = newdecay;
		Ptright     = old.Ptright - momentum->vect();
		Ptright.setZ(0.);
		LeftParton  = old.LeftParton;
		Ptleft      = old.Ptleft;
	} else
	{
		G4Exception("G4FragmentingString::G4FragmentingString: no decay Direction defined");
	}
	Pplus  = old.Pplus  - (momentum->e() + momentum->pz());
	Pminus = old.Pminus - (momentum->e() - momentum->pz());
	
	G4double Eold=0.5 * (old.Pplus + old.Pminus);
	G4double Enew=0.5 * (Pplus + Pminus);
}


//---------------------------------------------------------------------------------

G4FragmentingString::~G4FragmentingString()
{}


//---------------------------------------------------------------------------------

void G4FragmentingString::SetLeftPartonStable()
{
     theStableParton=GetLeftParton();
     theDecayParton=GetRightParton();
     decaying=Right;
}

//---------------------------------------------------------------------------------

void G4FragmentingString::SetRightPartonStable()
{
     theStableParton=GetRightParton();
     theDecayParton=GetLeftParton();
     decaying=Left;
}

//---------------------------------------------------------------------------------

G4int G4FragmentingString::GetDecayDirection() const
{
	if      (decaying == Left ) return +1;
	else if (decaying == Right) return -1;
	else G4Exception("G4FragmentingString::GetDecayDirection: decay side UNdefined!");
	return 0;
}
 
//---------------------------------------------------------------------------------

G4bool G4FragmentingString::FourQuarkString() const
{
	return   LeftParton->GetParticleSubType()== "di_quark" 
	     && RightParton->GetParticleSubType()== "di_quark";
}

//---------------------------------------------------------------------------------

G4bool G4FragmentingString::DecayIsQuark()
{
	return theDecayParton->GetParticleSubType()== "quark";
}

G4bool G4FragmentingString::StableIsQuark()
{
	return theStableParton->GetParticleSubType()== "quark";
}

//---------------------------------------------------------------------------------

G4ThreeVector G4FragmentingString::StablePt()
{
	if (decaying == Left ) return Ptright;
	else if (decaying == Right ) return Ptleft;
	else G4Exception("G4FragmentingString::DecayPt: decay side UNdefined!");
	return G4ThreeVector();
}

G4ThreeVector G4FragmentingString::DecayPt()
{
	if (decaying == Left ) return Ptleft;
	else if (decaying == Right ) return Ptright;
	else G4Exception("G4FragmentingString::DecayPt: decay side UNdefined!");
	return G4ThreeVector();
}

//---------------------------------------------------------------------------------

G4double G4FragmentingString::LightConePlus()
{
	return Pplus;
}

G4double G4FragmentingString::LightConeMinus()
{
	return Pminus;
}

G4double G4FragmentingString::LightConeDecay()
{
	if (decaying == Left ) return Pplus;
	else if (decaying == Right ) return Pminus;
	else G4Exception("G4FragmentingString::DecayPt: decay side UNdefined!");
	return 0;
}

//---------------------------------------------------------------------------------

G4LorentzVector G4FragmentingString::Get4Momentum() const
{
        G4LorentzVector momentum(Ptleft+Ptright,0);
	momentum.setPz(0.5*(Pplus-Pminus));
	momentum.setE(0.5*(Pplus+Pminus));
	return momentum;
}

G4double G4FragmentingString::Mass2() const
{
	return Pplus*Pminus - (Ptleft+Ptright).mag2();
}

G4double G4FragmentingString::Mass() const
{
	return sqrt(this->Mass2());
}
