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
#include "G4InelasticSplitableHadron.hh"

#include "G4ParticleDefinition.hh"
#include "G4Parton.hh"
#include "Randomize.hh"

G4InelasticSplitableHadron::G4InelasticSplitableHadron()
{
    udspin1 = 1./6.;
    uuspin1 = 1./3.;
    udspin0 = 1./2.;
}

G4InelasticSplitableHadron::G4InelasticSplitableHadron(const G4ReactionProduct & aPrimary)
      :  G4VSplitableHadron(aPrimary)
{
	PartonIndex=-2;
	Parton[0]=NULL;
}

G4InelasticSplitableHadron::G4InelasticSplitableHadron(const G4Nucleon & aNucleon)
      :  G4VSplitableHadron(aNucleon)
{
	PartonIndex=-2;
	Parton[0]=NULL;
}

G4InelasticSplitableHadron::~G4InelasticSplitableHadron()
{}

const G4InelasticSplitableHadron & G4InelasticSplitableHadron::operator=(const G4InelasticSplitableHadron &right)
{
  G4Exception("G4InelasticSplitableHadron::operator= meant to not be accessable");
  return *this;
}


//***********************************************************************************************

void G4InelasticSplitableHadron::SplitUp()
{
// Split once only...
	if (Parton[0] != NULL) return;

// flavours of quark ends
	G4int PDGcode=GetDefinition()->GetPDGEncoding();
	G4int stringStart, stringEnd;
	GetValenceQuarkFlavors(PDGcode, stringStart,stringEnd);
	Parton[0] = new G4Parton(stringStart);
	Parton[1] = new G4Parton(stringEnd);
	PartonIndex=-1;
	
}

G4Parton * G4InelasticSplitableHadron::GetNextParton()
{
	++PartonIndex;
	if ( PartonIndex > 1 || PartonIndex < 0 ) return NULL;
	
	return Parton[PartonIndex];
}

//
//----------------------- Implementation--------------------------
//
void G4InelasticSplitableHadron::GetValenceQuarkFlavors(G4int PDGcode, G4int& aEnd, G4int& bEnd)
   {
   // Note! convention aEnd = q or qqbar and bEnd = qbar or qq.
   G4int absPDGcode = abs(PDGcode);
   if (absPDGcode < 1000)   
      {
      G4int heavy =  absPDGcode/100;
      G4int light = (absPDGcode%100)/10;
      G4int anti  = 1 - 2*(G4std::max(heavy, light)%2);
      if (PDGcode < 0 ) anti = -anti;
      heavy *= anti;
      light *= -1 * anti;
      if ( heavy > 0) 
	 {
	 aEnd=heavy;
	 bEnd=light;
	 } 
      else 
	 {
	 aEnd=light;
	 bEnd=heavy;
	 }
      return; 	 
      }
   G4int j1000 = PDGcode/ 1000;
   G4int j100  = (PDGcode % 1000)/100;
   G4int j10   = (PDGcode % 100)/10;
   G4double random = G4UniformRand();
   if (abs(j100) >= abs(j10) )
      {  	    
      if ( random < udspin1 )
	 {
	 aEnd = j1000;
	 bEnd = j100 + j10 + 1;
	 } 
      else if ( random < (udspin1 + uuspin1) )
	 {
	 aEnd = j10;
	 bEnd = j1000 + j100 + 1;
	 } 
      else
	 {
	 aEnd = j100;
	 bEnd = j1000 + j10 + 0;
	 }
      if (aEnd < 0)
	 {
	 G4int Swap = aEnd;
	 aEnd = bEnd;
	 bEnd = Swap;
	 }
      return;
      } 
   // Lambda-like hadrons have two lightest quarks in spin 0
   if ( random < udspin1 )
      {
      aEnd = j1000;
      bEnd = j100 + j10 + 0;
      } 
   else if ( random < (udspin1 + uuspin1) )
      {
      aEnd = j10;
      bEnd = j1000 + j100 + 1;
      } 
   else
      {
      aEnd = j100;
      bEnd = j1000 + j10 + 1;
      }		
   if (aEnd < 0)
      {
      G4int Swap = aEnd;
      aEnd = bEnd;
      bEnd = Swap;
      }
   }   
	
//***********************************************************************************************

