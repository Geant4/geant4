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
	theStableParton=old.theStableParton;
	theDecayParton=old.theDecayParton;
	decaying=old.decaying;
        Pstring=old.Pstring;
        Pleft  =old.Pleft;
        Pright =old.Pright;
}

G4FragmentingString & G4FragmentingString::operator =(const G4FragmentingString &old)
{
   if (this != &old)
   {
     LeftParton=old.LeftParton;
     RightParton=old.RightParton;
     Ptleft=old.Ptleft;
     Ptright=old.Ptright;
     Pplus=old.Pplus;
     Pminus=old.Pminus;
     theStableParton=old.theStableParton;
     theDecayParton=old.theDecayParton;
     decaying=old.decaying;
     Pstring=old.Pstring;
     Pleft  =old.Pleft;
     Pright =old.Pright;
   }
   return *this;
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
	theStableParton=0;
	theDecayParton=0;

        if (excited.GetDirection() > 0) {decaying=Left; }
        else                            {decaying=Right;}

        Pleft  = excited.GetLeftParton()->Get4Momentum();
        Pright = excited.GetRightParton()->Get4Momentum();
        Pstring= Pleft + Pright;

        Pplus  = Pstring.plus();
        Pminus = Pstring.minus();	
}

//---------------------------------------------------------------------------------

G4FragmentingString::G4FragmentingString(const G4FragmentingString &old,
					 G4ParticleDefinition * newdecay,
					 const G4LorentzVector *momentum)
{
	decaying=None;
	// Momentum of produced hadron
        G4LorentzVector Momentum = G4LorentzVector(momentum->vect(),momentum->e());
             
        if ( old.decaying == Left )
        {
                RightParton= old.RightParton;
                Ptright    = old.Ptright;
                Pright     = old.Pright;

                LeftParton = newdecay;
                Ptleft     = old.Ptleft - momentum->vect();
                Ptleft.setZ(0.);
                Pleft = old.Pleft - Momentum;

                Pstring = Pleft + Pright;
                Pplus   = Pstring.plus();
                Pminus  = Pstring.minus();

                theDecayParton=GetLeftParton();
                theStableParton=GetRightParton();
                decaying = Left;
        } else if ( old.decaying == Right )
        {
	        RightParton = newdecay;
                Ptright     = old.Ptright - momentum->vect();
                Ptright.setZ(0.);
                Pright = old.Pright - Momentum;

                LeftParton  = old.LeftParton;
                Ptleft      = old.Ptleft;
                Pleft  = old.Pleft;

                Pstring = Pleft + Pright;
                Pplus   = Pstring.plus();
                Pminus  = Pstring.minus();

                theDecayParton=GetRightParton();
                theStableParton=GetLeftParton();
                decaying = Right;
        } else
	{
		throw G4HadronicException(__FILE__, __LINE__, 
                  "G4FragmentingString::G4FragmentingString: no decay Direction defined");
	}
}


//---------------------------------------------------------------------------------

G4FragmentingString::G4FragmentingString(const G4FragmentingString &old,  
					 G4ParticleDefinition * newdecay) 
{                                                                         
	decaying=None;                                                    

        Ptleft.setX(0.);  Ptleft.setY(0.);  Ptleft.setZ(0.);
        Ptright.setX(0.); Ptright.setY(0.); Ptright.setZ(0.);
        Pplus=0.; Pminus=0.;                               
        theStableParton=0; theDecayParton=0;              

        Pstring=G4LorentzVector(0.,0.,0.,0.);
        Pleft  =G4LorentzVector(0.,0.,0.,0.);
        Pright =G4LorentzVector(0.,0.,0.,0.);

	if ( old.decaying == Left )                                       
	{                                                                 
		RightParton= old.RightParton;                             
		LeftParton = newdecay;                                    
                decaying=Left;
	} else if ( old.decaying == Right )                               
	{                                                                 
		RightParton = newdecay;                                   
		LeftParton  = old.LeftParton;                             
                decaying=Right;
	} else                                                            
	{
		throw G4HadronicException(__FILE__, __LINE__, 
                  "G4FragmentingString::G4FragmentingString: no decay Direction defined");
	}
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
	else throw G4HadronicException(__FILE__, __LINE__, 
               "G4FragmentingString::GetDecayDirection: decay side UNdefined!");
	return 0;
}
 
//---------------------------------------------------------------------------------

G4bool G4FragmentingString::IsAFourQuarkString() const
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
	else throw G4HadronicException(__FILE__, __LINE__, "G4FragmentingString::DecayPt: decay side UNdefined!");
	return G4ThreeVector();
}

G4ThreeVector G4FragmentingString::DecayPt()
{
	if (decaying == Left ) return Ptleft;
	else if (decaying == Right ) return Ptright;
	else throw G4HadronicException(__FILE__, __LINE__, "G4FragmentingString::DecayPt: decay side UNdefined!");
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
	else throw G4HadronicException(__FILE__, __LINE__, "G4FragmentingString::DecayPt: decay side UNdefined!");
}

//---------------------------------------------------------------------------------

G4LorentzVector G4FragmentingString::Get4Momentum() const
{
        return Pstring;
}

G4double G4FragmentingString::Mass2() const
{
  return Pstring.mag2();
}

G4double G4FragmentingString::Mass() const
{
  return Pstring.mag();
}

G4double G4FragmentingString::MassT2() const
{
  return Pplus*Pminus;
}

G4LorentzVector G4FragmentingString::GetPstring()
{ 
  return Pstring;
}

G4LorentzVector G4FragmentingString::GetPleft()
{ 
  return Pleft;
}

G4LorentzVector G4FragmentingString::GetPright()
{ 
  return Pright;
}

G4LorentzRotation G4FragmentingString::TransformToAlignedCms()
{
     G4LorentzVector momentum = Pstring;
     G4LorentzRotation toAlignedCms(-1*momentum.boostVector());
     momentum = toAlignedCms * Pleft;

     toAlignedCms.rotateZ(-1*momentum.phi());
     toAlignedCms.rotateY(-1*momentum.theta());

     Pleft   *= toAlignedCms;
     Pright  *= toAlignedCms;
     Pstring *= toAlignedCms;

     Ptleft  = G4ThreeVector(Pleft.vect());
     Ptleft.setZ(0.);
     Ptright = G4ThreeVector(Pright.vect());
     Pplus  = Pstring.plus();
     Pminus = Pstring.minus();

     return toAlignedCms;
}
