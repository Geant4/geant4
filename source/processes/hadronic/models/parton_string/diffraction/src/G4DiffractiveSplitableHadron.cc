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
// $Id$
// GEANT4 tag $Name:  $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4DiffractiveSplitableHadron----------------
//             by Gunter Folger, August 1998.
//       class splitting an interacting particle. Used by FTF String Model.
// ------------------------------------------------------------

#include "G4DiffractiveSplitableHadron.hh"

#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

G4DiffractiveSplitableHadron::G4DiffractiveSplitableHadron()

{
	PartonIndex=-1;
	Parton[0] = new G4Parton(1);
	Parton[1] = new G4Parton(-1);
}

G4DiffractiveSplitableHadron::G4DiffractiveSplitableHadron(const G4ReactionProduct & aPrimary)
      :  G4VSplitableHadron(aPrimary)
{
	PartonIndex=-2;
	Parton[0]=NULL;
}

G4DiffractiveSplitableHadron::G4DiffractiveSplitableHadron(const G4Nucleon & aNucleon)
      :  G4VSplitableHadron(aNucleon)
{
	PartonIndex=-2;
	Parton[0]=NULL;
}

G4DiffractiveSplitableHadron::G4DiffractiveSplitableHadron(const G4VKineticNucleon * aNucleon)
      :  G4VSplitableHadron(aNucleon)
{
	PartonIndex=-2;
	Parton[0]=NULL;
}

G4DiffractiveSplitableHadron::~G4DiffractiveSplitableHadron()
{
//G4cout<<"Destruct G4DiffractiveSplitableHadron"<<Parton[0]<<" "<<Parton[1]<<G4endl;
//        if(Parton[0] != NULL){delete Parton[0]; delete Parton[1];}
}


void G4DiffractiveSplitableHadron::SplitUp()
{
//G4cout<<"SplitUp() IsSplit() Parton[0] "<<IsSplit()<<" "<<Parton[0]<<G4endl;
  if (IsSplit()) return;
  Splitting();
// Split once only...
	if (Parton[0] != NULL) return;

// flavours of quark ends
	
	G4int PDGcode=GetDefinition()->GetPDGEncoding();

	G4int stringStart, stringEnd;
	ChooseStringEnds(PDGcode, &stringStart,&stringEnd);

	Parton[0] = new G4Parton(stringStart);
	Parton[1] = new G4Parton(stringEnd);
	PartonIndex=-1;
}

G4Parton * G4DiffractiveSplitableHadron::GetNextParton()
{
	++PartonIndex;
	if ( PartonIndex > 1 || PartonIndex < 0 ) return NULL;
        G4int PartonInd(PartonIndex);        // Vova
        if(PartonIndex == 1) PartonIndex=-1; // Vova	
	return Parton[PartonInd];
//	return Parton[PartonIndex];
}

G4Parton * G4DiffractiveSplitableHadron::GetNextAntiParton()
{
	++PartonIndex;   // Uzhi 22.11.10
	if ( PartonIndex > 1 || PartonIndex < 0 ) return NULL;
        G4int PartonInd(PartonIndex);
        if(PartonIndex == 1) PartonIndex=-1;	
	return Parton[PartonInd];
//  return NULL; // to be looked at @@
}

void G4DiffractiveSplitableHadron::SetFirstParton(G4int PDGcode)// Uzhi 24.11.10
{
         delete Parton[0];
         Parton[0]=new G4Parton(PDGcode);
}

void G4DiffractiveSplitableHadron::SetSecondParton(G4int PDGcode)// Uzhi 24.11.10
{
         delete Parton[1];
         Parton[1]=new G4Parton(PDGcode);
}
//
//----------------------- Implementation--------------------------
//
void G4DiffractiveSplitableHadron::ChooseStringEnds(G4int PDGcode,G4int * aEnd, G4int * bEnd) const
{
	const G4double udspin1= 1./6.;
	const G4double uuspin1= 1./3.;
//	const G4double udspin0= 1./2.; //@
	
	G4int absPDGcode=std::abs(PDGcode);
	
	if ( absPDGcode < 1000 )   //--------------------  Meson -------------
	{
	   G4int heavy=  absPDGcode/ 100;
	   G4int light= (absPDGcode %100)/10;
	   
//	    G4int anti= std::pow(-1 , std::max( heavy, light));
	    G4int anti= 1 -2 * ( std::max( heavy, light ) % 2 );
	    if (PDGcode < 0 ) anti *=-1;
	    
	    heavy *= anti;
	    light *= -1 * anti;
	    
	    if ( G4UniformRand() < 0.5 ) 
	    {
	        *aEnd=heavy;
	        *bEnd=light;
	    } 
	    else 
	    {
	        *aEnd=light;
	        *bEnd=heavy;
	    }
	}
	else                      //-------------------- Barion --------------
	{
	    G4int j1000 = PDGcode/ 1000;
	    G4int j100  = (PDGcode % 1000) / 100;
	    G4int j10   = (PDGcode % 100) / 10;
	
	    G4double random= G4UniformRand();

	
	     if ( std::abs(j100) >= std::abs(j10) )
	     {  	    
	        if ( random < udspin1 )
	        {
	            *aEnd=j1000;
	            *bEnd= Diquark( j100, j10, 1);
	        } else if ( random < (udspin1 + uuspin1) )
	        {
	            *aEnd= j10;
	            *bEnd= Diquark( j1000, j100, 1);
	        } else
	        {
	            *aEnd=j100;
		       // Careful, there is no diquark q1=q2, (q1 q2)0, 
		       //   possible for Omega-
	            *bEnd= Diquark( j1000, j10, j1000 != j100 ? 0 : 1);
	        }
	     } else
	     {
// Lambda-like hadrons have two lightest quarks in spin 0
	        if ( random < udspin1 )
	        {
	            *aEnd=j1000;
		    		// as above, but with charmed barions
	            *bEnd= Diquark( j100, j10, j100 != j10 ? 0 : 10);
	        } else if ( random < (udspin1 + uuspin1) )
	        {
	            *aEnd= j10;
	            *bEnd= Diquark( j1000, j100, 1);
	        } else
	        {
	            *aEnd=j100;
	            *bEnd= Diquark( j1000, j10, 1);
	        }
		
	     } 
	
	}   
	    
	     
	    
}

G4int G4DiffractiveSplitableHadron::Diquark(G4int aquark,G4int bquark,G4int Spin) const
{
	G4int diquarkPDG;
	
	diquarkPDG = std::max( std::abs(aquark), std::abs(bquark) )*1000 +
	             std::min( std::abs(aquark), std::abs(bquark) )*100  +
	             2*Spin + 1;
	return ( aquark > 0 && bquark > 0 ) ?  diquarkPDG : -1*diquarkPDG;
}	      
