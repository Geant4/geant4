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
// $Id: G4QGSDiffractiveExcitation.cc 102316 2017-01-20 16:12:52Z gcosmo $
// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      ---------------- G4QGSDiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//         diffractive Excitation used by strings models
//	   Take a projectile and a target
//	   excite the projectile and target
//  Essential changed by V. Uzhinsky in November - December 2006
//  in order to put it in a correspondence with original FRITIOF
//  model. Variant of FRITIOF with nucleon de-excitation is implemented.  
// ---------------------------------------------------------------------

// Modified:
//  25-05-07 : G.Folger
//       move from management/G4DiffractiveExcitation to to qgsm/G4QGSDiffractiveExcitation
//                  

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4QGSDiffractiveExcitation.hh"
#include "G4LorentzRotation.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"
//#include "G4ios.hh"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

//============================================================================

//#define debugDoubleDiffraction

//============================================================================

G4QGSDiffractiveExcitation::G4QGSDiffractiveExcitation()
{
}

G4QGSDiffractiveExcitation::~G4QGSDiffractiveExcitation()
{
}


G4bool G4QGSDiffractiveExcitation::
ExciteParticipants(G4VSplitableHadron *projectile, G4VSplitableHadron *target, G4bool ) const  // Uzhi Oct. 2016 , G4bool ProjectileDiffraction
{
  #ifdef debugDoubleDiffraction
     G4cout<<G4endl<<"G4QGSDiffractiveExcitation::ExciteParticipants - Double diffraction."<<G4endl;
     G4cout<<"Proj Targ "<<projectile->GetDefinition()->GetParticleName()<<" "<<target->GetDefinition()->GetParticleName()<<G4endl;
     G4cout<<"Proj 4 Mom "<<projectile->Get4Momentum()<<" "<<projectile->Get4Momentum().mag()<<G4endl;
     G4cout<<"Targ 4 Mom "<<target->Get4Momentum()    <<" "<<target->Get4Momentum().mag()    <<G4endl;
  #endif



	G4LorentzVector Pprojectile=projectile->Get4Momentum();

	// -------------------- Projectile parameters -----------------------------------
	G4bool PutOnMassShell=0;

	G4double M0projectile = Pprojectile.mag();                        // Without de-excitation

	if(M0projectile < projectile->GetDefinition()->GetPDGMass())
	{
		PutOnMassShell=1;
		M0projectile=projectile->GetDefinition()->GetPDGMass();
	}

	// -------------------- Target parameters ----------------------------------------------
	G4LorentzVector Ptarget=target->Get4Momentum();

	G4double M0target = Ptarget.mag();

	if(M0target < target->GetDefinition()->GetPDGMass())
	{
		PutOnMassShell=1;
		M0target=target->GetDefinition()->GetPDGMass();
	}

	G4LorentzVector Psum=Pprojectile+Ptarget;
	G4double S=Psum.mag2();
	G4double SqrtS=std::sqrt(S);

	if(SqrtS < M0projectile + M0target) {return false;} // The model cannot work for pp-interactions
	                                                    // at Plab < 1.3 GeV/c. Uzhi

	G4double Mprojectile2 = M0projectile * M0projectile;
	G4double Mtarget2     = M0target     * M0target;    //Ptarget.mag2(); // for AA-inter.

	// Transform momenta to cms and then rotate parallel to z axis;

	G4LorentzRotation toCms(-1*Psum.boostVector());

	G4LorentzVector Ptmp=toCms*Pprojectile;

	if ( Ptmp.pz() <= 0. )
	{
		// "String" moving backwards in  CMS, abort collision !!
		//G4cout << " abort Collision!! " << G4endl;
		return false;
	}

	toCms.rotateZ(-1*Ptmp.phi());
	toCms.rotateY(-1*Ptmp.theta());

	G4LorentzRotation toLab(toCms.inverse());

	Pprojectile.transform(toCms);
	Ptarget.transform(toCms);

	G4double PZcms2=(S*S+Mprojectile2*Mprojectile2+Mtarget2*Mtarget2-
			 2*S*Mprojectile2-2*S*Mtarget2-2*Mprojectile2*Mtarget2)/4./S;

	if(PZcms2 < 0) {return false;}   // It can be in an interaction with off-shell nuclear nucleon

	G4double PZcms = std::sqrt(PZcms2);

	if(PutOnMassShell)
	{
		if(Pprojectile.z() > 0.)
		{
			Pprojectile.setPz( PZcms);
			Ptarget.setPz(    -PZcms);
		}
		else
		{
			Pprojectile.setPz(-PZcms);
			Ptarget.setPz(     PZcms);
		};

		Pprojectile.setE(std::sqrt(Mprojectile2+sqr(Pprojectile.x())+sqr(Pprojectile.y())+PZcms2));
		Ptarget.setE(    std::sqrt(    Mtarget2+sqr(    Ptarget.x())+sqr(    Ptarget.y())+PZcms2));
	}

	G4double maxPtSquare = PZcms2;

  #ifdef debugDoubleDiffraction
	G4cout << "Pprojectile after boost to CMS: " << Pprojectile <<" "<<Pprojectile.mag()<<G4endl;
	G4cout << "Ptarget     after boost to CMS: " << Ptarget     <<" "<<Ptarget.mag()    <<G4endl;
  #endif

	G4int    PrPDGcode=projectile->GetDefinition()->GetPDGEncoding();
	G4int    absPrPDGcode=std::abs(PrPDGcode);
	G4double MinPrDiffMass(0.);
	G4double AveragePt2(0.);

	if(M0projectile <= projectile->GetDefinition()->GetPDGMass())
        { // Normal projectile 
	  if( absPrPDGcode > 1000 )                         //------Projectile is baryon --------
	  {
		MinPrDiffMass = 1.16;                       // GeV
		AveragePt2 = 0.3;                           // GeV^2
	  }
	  else if( absPrPDGcode == 211 || PrPDGcode ==  111) //------Projectile is Pion -----------
	  {
		MinPrDiffMass = 1.0;                       // GeV
		AveragePt2 = 0.3;                          // GeV^2
	  }
	  else if( absPrPDGcode == 321 || absPrPDGcode == 130 || absPrPDGcode == 310) //-Projectile is Kaon-
	  {
		MinPrDiffMass = 1.1;                        // GeV
		AveragePt2 = 0.3;                           // GeV^2
	  }
	  else                                              //------Projectile is undefined, Nucleon assumed
	  {
		MinPrDiffMass = 1.16;                       // GeV
		AveragePt2 = 0.3;                           // GeV^2
	  }
        }
        else
        { // Excited projectile
		MinPrDiffMass = M0projectile + 220.0*MeV;
		AveragePt2 = 0.3; 
        }

	MinPrDiffMass = MinPrDiffMass * GeV;
	AveragePt2 = AveragePt2 * GeV*GeV;
//---------------------------------------------
	G4double MinTrDiffMass = 1.16*GeV;

	if(SqrtS < MinPrDiffMass + MinTrDiffMass) {return false;}   // The model cannot work at low energy

	G4double MinPrDiffMass2 = MinPrDiffMass * MinPrDiffMass;
	G4double MinTrDiffMass2 = MinTrDiffMass * MinTrDiffMass;

	G4double Pt2;
	G4double ProjMassT2, ProjMassT;
	G4double TargMassT2, TargMassT;
	G4double PMinusNew, TPlusNew;

	G4LorentzVector Qmomentum;
	G4double Qminus, Qplus;

	G4int whilecount=0;
	do {
		if (whilecount++ >= 500 && (whilecount%100)==0)
			//	   	 G4cout << "G4QGSDiffractiveExcitation::ExciteParticipants possibly looping"
			//	   	 << ", loop count/ maxPtSquare : "
			//           	 << whilecount << " / " << maxPtSquare << G4endl;
			if (whilecount > 1000 )
			{
				Qmomentum=G4LorentzVector(0.,0.,0.,0.);
				return false; 	  //  Ignore this interaction
			}

		//  Generate pt
		Qmomentum=G4LorentzVector(GaussianPt(AveragePt2,maxPtSquare),0);

		Pt2=G4ThreeVector(Qmomentum.vect()).mag2();
		ProjMassT2=MinPrDiffMass2+Pt2;
		ProjMassT =std::sqrt(ProjMassT2);

		TargMassT2=MinTrDiffMass2+Pt2;
		TargMassT =std::sqrt(TargMassT2);

                if(SqrtS < ProjMassT + TargMassT) continue;
 
		PZcms2=(S*S+ProjMassT2*ProjMassT2+
				TargMassT2*TargMassT2-
				2.*S*ProjMassT2-2.*S*TargMassT2-
				2.*ProjMassT2*TargMassT2)/4./S;
		if(PZcms2 < 0 ) {PZcms2=0;};
		PZcms =std::sqrt(PZcms2);

		G4double PMinusMin=std::sqrt(ProjMassT2+PZcms2)-PZcms;
		G4double PMinusMax=SqrtS-TargMassT;

		PMinusNew=ChooseP(PMinusMin,PMinusMax);
		Qminus=PMinusNew-Pprojectile.minus();

		G4double TPlusMin=std::sqrt(TargMassT2+PZcms2)-PZcms;
		G4double TPlusMax=SqrtS-ProjMassT;

		TPlusNew=ChooseP(TPlusMin, TPlusMax);
		Qplus=-(TPlusNew-Ptarget.plus());

		Qmomentum.setPz( (Qplus-Qminus)/2 );
		Qmomentum.setE(  (Qplus+Qminus)/2 );

	} while ( (Pprojectile+Qmomentum).mag2() <  MinPrDiffMass2 ||      // Uzhi No without excitation
		  (Ptarget    -Qmomentum).mag2() <  MinTrDiffMass2   ); 

	Pprojectile += Qmomentum;
	Ptarget     -= Qmomentum;

	// Transform back and update SplitableHadron Participant.
	Pprojectile.transform(toLab);
	Ptarget.transform(toLab);

  #ifdef debugDoubleDiffraction
	G4cout << "Pprojectile after boost to Lab: " << Pprojectile <<" "<<Pprojectile.mag()<<G4endl;
	G4cout << "Ptarget     after boost to Lab: " << Ptarget     <<" "<<Ptarget.mag()    <<G4endl;
  #endif

	target->Set4Momentum(Ptarget);
	projectile->Set4Momentum(Pprojectile);

	return true;
}


G4ExcitedString * G4QGSDiffractiveExcitation::
String(G4VSplitableHadron * hadron, G4bool isProjectile) const
{
	hadron->SplitUp();
	G4Parton *start= hadron->GetNextParton();
	if ( start==NULL) 
	{ G4cout << " G4QGSDiffractiveExcitation::String() Error:No start parton found"<< G4endl;
	return NULL;
	}
	G4Parton *end  = hadron->GetNextParton();
	if ( end==NULL) 
	{ G4cout << " G4QGSDiffractiveExcitation::String() Error:No end parton found"<< G4endl;
	return NULL;
	}

	G4ExcitedString * string;
	if ( isProjectile ) 
	{
		string= new G4ExcitedString(end,start, +1);
	} else {
		string= new G4ExcitedString(start,end, -1);
	}

	string->SetPosition(hadron->GetPosition());

	// momenta of string ends
/*                                                                  // Uzhi 2016
	G4double ptSquared= hadron->Get4Momentum().perp2();
	G4double transverseMassSquared= hadron->Get4Momentum().plus()
				    		*	hadron->Get4Momentum().minus();


	G4double maxAvailMomentumSquared=
			sqr( std::sqrt(transverseMassSquared) - std::sqrt(ptSquared) );
*/
	G4double maxAvailMomentumSquared=sqr(hadron->Get4Momentum().mag()/2.);  // Uzhi 2016

	G4double widthOfPtSquare = 0.5*sqr(GeV);   //0.25; // Uzhi 2016 // Uzhi <Pt^2>=0.25 ??????????????????
	G4ThreeVector pt=GaussianPt(widthOfPtSquare,maxAvailMomentumSquared);

	G4LorentzVector Pstart(G4LorentzVector(pt,0.));
	G4LorentzVector Pend;
	Pend.setPx(hadron->Get4Momentum().px() - pt.x());
	Pend.setPy(hadron->Get4Momentum().py() - pt.y());

	G4double tm1=hadron->Get4Momentum().minus() +
			( Pend.perp2()-Pstart.perp2() ) / hadron->Get4Momentum().plus();

	G4double tm2= std::sqrt( std::max(0., sqr(tm1) -
			4. * Pend.perp2() * hadron->Get4Momentum().minus()
			/  hadron->Get4Momentum().plus() ));

	G4int Sign= isProjectile ? -1 : 1;

	G4double endMinus  = 0.5 * (tm1 + Sign*tm2);
	G4double startMinus= hadron->Get4Momentum().minus() - endMinus;

	G4double startPlus= Pstart.perp2() /  startMinus;
	G4double endPlus  = hadron->Get4Momentum().plus() - startPlus;

	Pstart.setPz(0.5*(startPlus - startMinus));
	Pstart.setE(0.5*(startPlus + startMinus));

	Pend.setPz(0.5*(endPlus - endMinus));
	Pend.setE(0.5*(endPlus + endMinus));

	start->Set4Momentum(Pstart);
	end->Set4Momentum(Pend);

	#ifdef debugQGSdiffExictation
		G4cout << " generated string flavors          " << start->GetPDGcode()   << " / " << end->GetPDGcode() << G4endl;
		G4cout << " generated string momenta:   quark " << start->Get4Momentum() << "mass : " <<start->Get4Momentum().mag()<< G4endl;
		G4cout << " generated string momenta: Diquark " << end ->Get4Momentum()  << "mass : " <<end->Get4Momentum().mag()<< G4endl;
		G4cout << " sum of ends                       " << Pstart+Pend << G4endl;
		G4cout << " Original                          " << hadron->Get4Momentum() << G4endl;
	#endif

	return string;	
}


// --------- private methods ----------------------

G4double G4QGSDiffractiveExcitation::ChooseP(G4double Pmin, G4double Pmax) const
{
	// choose an x between Xmin and Xmax with P(x) ~ 1/x
	//  to be improved...

	G4double range=Pmax-Pmin;

	if ( Pmin <= 0. || range <=0. ) 
	{
		G4cout << " Pmin, range : " << Pmin << " , " << range << G4endl;
		throw G4HadronicException(__FILE__, __LINE__, "G4QGSDiffractiveExcitation::ChooseP : Invalid arguments ");
	}

	G4double P;
	P=Pmin * G4Pow::GetInstance()->powA(Pmax/Pmin,G4UniformRand());
	//debug-hpw	cout << "DiffractiveX "<<x<<G4endl;
	return P;
}

G4ThreeVector G4QGSDiffractiveExcitation::GaussianPt(G4double AveragePt2, G4double maxPtSquare) const
{            //  @@ this method is used in FTFModel as well. Should go somewhere common!

	G4double Pt2;

	Pt2 = -AveragePt2 * G4Log(1. + G4UniformRand() * (G4Exp(-maxPtSquare/AveragePt2)-1.));

	G4double Pt=std::sqrt(Pt2);

	G4double phi=G4UniformRand() * twopi;

	return G4ThreeVector (Pt*std::cos(phi), Pt*std::sin(phi), 0.);    
}
