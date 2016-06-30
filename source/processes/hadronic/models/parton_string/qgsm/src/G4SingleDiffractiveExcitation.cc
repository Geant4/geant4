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
// $Id: G4SingleDiffractiveExcitation.cc 94750 2015-12-07 08:24:29Z gcosmo $
// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      ---------------- G4SingleDiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// ------------------------------------------------------------

#include "G4SingleDiffractiveExcitation.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4LorentzRotation.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"

#include "G4Log.hh"
#include "G4Pow.hh"


G4SingleDiffractiveExcitation::G4SingleDiffractiveExcitation(G4double sigmaPt, G4double minextraMass,G4double x0mass)
:
widthOfPtSquare(-2*sqr(sigmaPt)) , minExtraMass(minextraMass),
minmass(x0mass)
{}

G4SingleDiffractiveExcitation::~G4SingleDiffractiveExcitation()
{}

G4bool G4SingleDiffractiveExcitation::
ExciteParticipants(G4VSplitableHadron *projectile, G4VSplitableHadron *target) const
{
/*
G4cout<<G4endl<<"G4SingleDiffractiveExcitation::ExciteParticipants"<<G4endl;
G4cout<<"Proj Targ "<<projectile->GetDefinition()->GetPDGEncoding()<<" "<<target->GetDefinition()->GetPDGEncoding()<<G4endl;
G4cout<<"minExtraMass "<<minExtraMass<<" minmass "<<minmass<<" widthOfPtSquare "<<widthOfPtSquare<<G4endl;
*/
	G4LorentzVector Pprojectile=projectile->Get4Momentum();
	G4double Mprojectile =    projectile->GetDefinition()->GetPDGMass();
	G4double Mprojectile2=sqr(projectile->GetDefinition()->GetPDGMass()); // + minExtraMass);

	G4LorentzVector Ptarget=target->Get4Momentum();
	G4double Mtarget =    target->GetDefinition()->GetPDGMass();
	G4double Mtarget2=sqr(target->GetDefinition()->GetPDGMass()); // + minExtraMass);

//G4cout<<"Pr Tr 4-Mom "<<Pprojectile<<" "<<Pprojectile.mag()<<G4endl<<"            "<<Ptarget    <<" "<<Ptarget.mag()   <<G4endl;

	// Transform momenta to cms and then rotate parallel to z axis;
G4double AveragePt2=sqr(400.*MeV);

	G4LorentzVector Psum=Pprojectile+Ptarget;
        G4double SqrtS=Psum.mag();
        G4double S    =Psum.mag2();

        if(SqrtS-Mprojectile-Mtarget <= 250.0*MeV) {
return true;
/*
           G4cerr<<"Projectile: "<<projectile->GetDefinition()->GetPDGEncoding()<<" "
                <<Pprojectile<<" "<<Pprojectile.mag()<<G4endl;
           G4cerr<<"Target:     "<<target->GetDefinition()->GetPDGEncoding()<<" "
                <<Ptarget<<" "<<Ptarget.mag()<<G4endl; 
           G4cerr<<"sqrt(S) = "<<SqrtS<<" Mp + Mt = "<<Pprojectile.mag()+Ptarget.mag()<<G4endl;
           throw G4HadronicException(__FILE__, __LINE__, "The QGSM cannot work at such low energy!"); 
*/
        }

	G4LorentzRotation toCms(-1*Psum.boostVector());

	G4LorentzVector Ptmp=toCms*Pprojectile;

	if ( Ptmp.pz() <= 0. )
	{
		// "String" moving backwards in  CMS, abort collision !!
		//	   	   G4cout << " abort Collision!! " << G4endl;
		return false;
	}

	toCms.rotateZ(-1*Ptmp.phi());
	toCms.rotateY(-1*Ptmp.theta());

	G4LorentzRotation toLab(toCms.inverse());

//G4cout << "Pprojectile  be4 boost " << Pprojectile << G4endl;
//G4cout << "Ptarget be4 boost :    " << Ptarget << G4endl;
	Pprojectile.transform(toCms);
	Ptarget.transform(toCms);
//G4cout << "Pprojectile  aft boost " << Pprojectile << G4endl;
//G4cout << "Ptarget aft boost :    " << Ptarget << G4endl;
        G4double maxPtSquare=sqr(Ptarget.pz());


	G4double Pt2, PZcms, PZcms2;
	G4double ProjMassT2, ProjMassT;
	G4double TargMassT2, TargMassT;
	G4double PMinusMin, PMinusMax;
	//G4double PPlusMin , PPlusMax;
	G4double TPlusMin, TPlusMax;
	G4double PMinusNew, PPlusNew, TPlusNew, TMinusNew;

	G4LorentzVector Qmomentum;
	G4double Qminus, Qplus;

	G4bool ProjectileDiffraction= G4UniformRand() > 0.5;
	if ( ProjectileDiffraction )
	{       // The projectile will fragment, the target will saved.
		Mprojectile2=sqr(Mprojectile   + 250.*MeV );

	} else {// The target will fragment, the projectile will saved.
		Mtarget2 = sqr(Mtarget + 250.*MeV );
	}

	G4int whilecount=0;
	do {
		whilecount++;

		if (whilecount > 1000 )
		{
//G4cout<<"whilecount > 1000 "<<whilecount<<G4endl;
			Qmomentum=G4LorentzVector(0.,0.,0.,0.);
			return false; 	  //  Ignore this interaction
		}
		
		//  Generate pt
		Qmomentum=G4LorentzVector(GaussianPt(AveragePt2,maxPtSquare),0);

		Pt2 = G4ThreeVector( Qmomentum.vect() ).mag2();
		ProjMassT2 = Mprojectile2 + Pt2;
		ProjMassT = std::sqrt( ProjMassT2 );
		TargMassT2 = Mtarget2 + Pt2;
		TargMassT = std::sqrt( TargMassT2 );
//G4cout<<whilecount<<" "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<SqrtS<<" "<<S<<" "<<ProjectileDiffraction<<G4endl;

		if ( SqrtS < ProjMassT + TargMassT ) continue;

		PZcms2 = ( S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2
                   - 2.0*S*ProjMassT2 - 2.0*S*TargMassT2 - 2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

		if ( PZcms2 < 0 ) continue;

		PZcms = std::sqrt( PZcms2 );

		if ( ProjectileDiffraction )
		{       // The projectile will fragment, the target will saved.
			PMinusMin = std::sqrt( ProjMassT2 + PZcms2 ) - PZcms;
        		PMinusMax = SqrtS - TargMassT;

			PMinusNew = ChooseX( PMinusMin, PMinusMax );
			TMinusNew = SqrtS - PMinusNew;

	        	Qminus = Ptarget.minus() - TMinusNew;
		        TPlusNew = TargMassT2 / TMinusNew;
	        	Qplus = Ptarget.plus() - TPlusNew;

		} else {// The target will fragment, the projectile will saved.
		        TPlusMin = std::sqrt( TargMassT2 + PZcms2 ) - PZcms;
		        TPlusMax = SqrtS - ProjMassT;
	
		        TPlusNew = ChooseX( TPlusMin, TPlusMax );
		        PPlusNew = SqrtS - TPlusNew;

		        Qplus = PPlusNew - Pprojectile.plus();
		        PMinusNew = ProjMassT2 / PPlusNew;
		        Qminus = PMinusNew - Pprojectile.minus();
		}


	        Qmomentum.setPz( (Qplus - Qminus)/2 );
	        Qmomentum.setE(  (Qplus + Qminus)/2 );
//G4cout<<ProjectileDiffraction<<" "<<( Pprojectile + Qmomentum ).mag2()<<" "<< Mprojectile2<<G4endl;
//G4cout<<!ProjectileDiffraction<<" "<<( Ptarget     + Qmomentum ).mag2()<<" "<< Mtarget2<<G4endl;
	} while ( ( ProjectileDiffraction&&( Pprojectile + Qmomentum ).mag2() <  Mprojectile2 ) ||
                  (!ProjectileDiffraction&&( Ptarget     - Qmomentum ).mag2() <  Mtarget2       )   );  /* Loop checking, 07.08.2015, A.Ribon */
                // Repeat the sampling because there was not any excitation


	Pprojectile += Qmomentum;

	Ptarget     -= Qmomentum;

	// Transform back and update SplitableHadron Participant.
	Pprojectile.transform(toLab);
	Ptarget.transform(toLab);

//G4cout << "Pprojectile  aft boost " << Pprojectile << G4endl;
//G4cout << "Ptarget aft boost :    " << Ptarget << G4endl;
//G4cout << "G4SingleDiffractiveExcitation- Target mass      " <<  Ptarget.mag() << G4endl;
//G4cout << "G4SingleDiffractiveExcitation- Projectile mass  " <<  Pprojectile.mag() << G4endl;

//G4int Uzhi; G4cin>>Uzhi;

	target->Set4Momentum(Ptarget);
	projectile->Set4Momentum(Pprojectile);

	return true;
}




// --------- private methods ----------------------

G4double G4SingleDiffractiveExcitation::ChooseX(G4double Xmin, G4double Xmax) const
{
	// choose an x between Xmin and Xmax with P(x) ~ 1/x
	G4double range=Xmax-Xmin;

	if ( Xmin <= 0. || range <=0. ) 
	{
		G4cout << " Xmin, range : " << Xmin << " , " << range << G4endl;
		throw G4HadronicException(__FILE__, __LINE__, "G4SingleDiffractiveExcitation::ChooseX : Invalid arguments ");
	}

	G4double x = Xmin*G4Pow::GetInstance()->powA(Xmax/Xmin, G4UniformRand() );
	return x;
}

G4ThreeVector G4SingleDiffractiveExcitation::GaussianPt(G4double widthSquare, G4double maxPtSquare) const
{            //  @@ this method is used in FTFModel as well. Should go somewhere common!

	G4double pt2;

        const G4int maxNumberOfLoops = 1000;
        G4int loopCounter = 0;
	do {
		pt2=-widthSquare * G4Log( G4UniformRand() );
	} while ( ( pt2 > maxPtSquare) && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
        if ( loopCounter >= maxNumberOfLoops ) {
          pt2 = 0.99*maxPtSquare;  // Just an acceptable value, without any physics consideration. 
        }

	pt2=std::sqrt(pt2);

	G4double phi=G4UniformRand() * twopi;

	return G4ThreeVector (pt2*std::cos(phi), pt2*std::sin(phi), 0.);    
}





