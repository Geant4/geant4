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
// $Id: G4QuarkExchange.cc 99348 2016-09-19 08:39:04Z vuzhinsk $
// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      ---------------- G4QuarkExchange --------------
//             by V. Uzhinsky, October 2016.
//       QuarkExchange is used by strings models.
//	    Take a projectile and a target.
//Simulate Q exchange with excitation of projectile or target.
// ------------------------------------------------------------

#include "G4QuarkExchange.hh"
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

//#define debugQuarkExchange

G4QuarkExchange::G4QuarkExchange(){}

G4QuarkExchange::~G4QuarkExchange(){}

G4bool G4QuarkExchange::
ExciteParticipants(G4VSplitableHadron *projectile, G4VSplitableHadron *target) const
{
  #ifdef debugQuarkExchange
     G4cout<<G4endl<<"G4QuarkExchange::ExciteParticipants"<<G4endl;
  #endif

	G4LorentzVector Pprojectile = projectile->Get4Momentum();
	G4double Mprojectile        = projectile->GetDefinition()->GetPDGMass();
	G4double Mprojectile2       = sqr(Mprojectile);

	G4LorentzVector Ptarget = target->Get4Momentum();
	G4double Mtarget        = target->GetDefinition()->GetPDGMass();
	G4double Mtarget2       = sqr(Mtarget);

  #ifdef debugQuarkExchange
     G4cout<<"Proj Targ "<<projectile->GetDefinition()->GetPDGEncoding()<<" "<<target->GetDefinition()->GetPDGEncoding()<<G4endl;
     G4cout<<"Proj. 4-Mom "<<Pprojectile<<" "<<Pprojectile.mag()<<G4endl
           <<"Targ. 4-Mom "<<Ptarget    <<" "<<Ptarget.mag()   <<G4endl;
  #endif

	G4LorentzVector Psum=Pprojectile+Ptarget;
        G4double SqrtS=Psum.mag();
        G4double S    =Psum.mag2();

  #ifdef debugQuarkExchange
     G4cout<<"SS Mpr Mtr SqrtS-Mprojectile-Mtarget "<<SqrtS<<" "<<Mprojectile<<" "<<Mtarget
           <<" "<<SqrtS-Mprojectile-Mtarget<<G4endl;
  #endif
        if(SqrtS-Mprojectile-Mtarget <= 250.0*MeV) {
  #ifdef debugQuarkExchange
           G4cerr<<"Energy is too small for quark exchange!"<<G4endl;
           G4cerr<<"Projectile: "<<projectile->GetDefinition()->GetPDGEncoding()<<" "
                <<Pprojectile<<" "<<Pprojectile.mag()<<G4endl;
           G4cerr<<"Target:     "<<target->GetDefinition()->GetPDGEncoding()<<" "
                <<Ptarget<<" "<<Ptarget.mag()<<G4endl; 
           G4cerr<<"sqrt(S) = "<<SqrtS<<" Mp + Mt = "<<Pprojectile.mag()+Ptarget.mag()<<G4endl;
  #endif
           return true;
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

	Pprojectile.transform(toCms);
	Ptarget.transform(toCms);

  #ifdef debugQuarkExchange
     G4cout << "Pprojectile  in CMS " << Pprojectile << G4endl;
     G4cout << "Ptarget      in CMS " << Ptarget     << G4endl;
  #endif
        G4double maxPtSquare=sqr(Ptarget.pz());

	G4double ProjectileMinDiffrMass(0.), TargetMinDiffrMass(0.);
	G4double AveragePt2(0.);

	G4int              PDGcode=projectile->GetDefinition()->GetPDGEncoding();
        G4int           absPDGcode=std::abs(PDGcode); 

	G4bool ProjectileDiffraction = true;

	if( absPDGcode > 1000 ) { ProjectileDiffraction = G4UniformRand() <= 0.5; }
	if( (absPDGcode == 211) || (absPDGcode == 111) ) { ProjectileDiffraction = G4UniformRand() <= 0.66; }
	if( (absPDGcode == 321) || (absPDGcode == 311)  || 
            (   PDGcode == 130) || (   PDGcode == 310) ) { ProjectileDiffraction = G4UniformRand() <= 0.5; } // Vova ???

	if ( ProjectileDiffraction ) {
	  if( absPDGcode > 1000 )                            //------Projectile is baryon --------
	  {
		ProjectileMinDiffrMass = 1.16;              // GeV
		AveragePt2 = 0.3;                           // GeV^2
	  }
	  else if( absPDGcode == 211 || absPDGcode ==  111) //------Projectile is Pion -----------
	  {
		ProjectileMinDiffrMass = 1.0;               // GeV
		AveragePt2 = 0.3;                           // GeV^2
	  }
	  else if( absPDGcode == 321 || absPDGcode == 130 || absPDGcode == 310) //Projectile is Kaon
	  {
		ProjectileMinDiffrMass = 1.1;               // GeV
		AveragePt2 = 0.3;                           // GeV^2
	  }
	  else if( absPDGcode == 22)                        //------Projectile is Gamma -----------
	  {
		ProjectileMinDiffrMass = 0.25;             // GeV
		AveragePt2 = 0.36;                         // GeV^2
	  }
	  else                                             //------Projectile is undefined, Nucleon assumed
	  {
		ProjectileMinDiffrMass = 1.1;              // GeV
		AveragePt2 = 0.3;                          // GeV^2
	  };

	  ProjectileMinDiffrMass = ProjectileMinDiffrMass * GeV;
          Mprojectile2=sqr(ProjectileMinDiffrMass);
        }
        else
        {
          TargetMinDiffrMass = 1.16*GeV;                     // For target nucleon
          Mtarget2 = sqr( TargetMinDiffrMass) ;
	  AveragePt2 = 0.3;                                  // GeV^2
        }   // end of if ( ProjectileDiffraction )

	AveragePt2 = AveragePt2 * GeV*GeV;                   // Uzhi 6 Oct. 2016

//----------------------- 
	G4double Pt2, PZcms, PZcms2;
	G4double ProjMassT2, ProjMassT;
	G4double TargMassT2, TargMassT;
	G4double PMinusMin, PMinusMax,  sqrtPMinusMin, sqrtPMinusMax;
      //G4double PPlusMin , PPlusMax;
	G4double TPlusMin, TPlusMax,    sqrtTPlusMin,  sqrtTPlusMax;
	G4double PMinusNew, PPlusNew, TPlusNew(0.), TMinusNew;

	G4LorentzVector Qmomentum;
	G4double Qminus, Qplus;

	G4double x(0.), y(0.);
	G4int whilecount=0;
	do {
		whilecount++;

		if (whilecount > 1000 )
		{
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

  #ifdef debugQuarkExchange
     G4cout<<"whilecount  Pt2  ProjMassT  TargMassT  SqrtS  S  ProjectileDiffraction"<<G4endl;
     G4cout<<whilecount<<" "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<SqrtS<<" "<<S<<" "<<ProjectileDiffraction<<G4endl;
  #endif
		if ( SqrtS < ProjMassT + TargMassT ) continue;

		PZcms2 = ( S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2
                   - 2.0*S*ProjMassT2 - 2.0*S*TargMassT2 - 2.0*ProjMassT2*TargMassT2 ) / 4.0 / S;

		if ( PZcms2 < 0 ) continue;

		PZcms = std::sqrt( PZcms2 );

		if ( ProjectileDiffraction )
		{// The projectile will fragment, the target will saved.
		 PMinusMin = std::sqrt( ProjMassT2 + PZcms2 ) - PZcms;
        	 PMinusMax = SqrtS - TargMassT;
		 sqrtPMinusMin = std::sqrt(PMinusMin); sqrtPMinusMax = std::sqrt(PMinusMax);
//===========
		 if( absPDGcode > 1000 ) 
		 { 
		  while(true)
		  {
		   x=PMinusMax-(PMinusMax-PMinusMin)*G4UniformRand();
		   y=G4UniformRand();
		   if( y < G4Pow::GetInstance()->powA(1.0-x/PMinusMax,3.) ) break;
		  }
		  PMinusNew = x;
		 } else if( (absPDGcode == 211) || (absPDGcode == 111) ) 
		 {
		  while(true)
		  {
		   x=sqrtPMinusMax-(sqrtPMinusMax-sqrtPMinusMin)*G4UniformRand();
		   y=G4UniformRand();
		   if( y < 1.0-0.7 * x/sqrtPMinusMax ) break;                  // 0.7 for Pi Vova
		  }
		  PMinusNew = sqr(x);
		 } else if( (absPDGcode == 321) || (absPDGcode == 311)  || 
			    (   PDGcode == 130) || (   PDGcode == 310) ) 
		 { // For K-mesons it must be found !!! Vova
		  while(true)
		  {
		   x=sqrtPMinusMax-(sqrtPMinusMax-sqrtPMinusMin)*G4UniformRand();
		   y=G4UniformRand();
		   if( y < 1.0-0.7 * x/sqrtPMinusMax ) break;
		  }
		  PMinusNew = sqr(x);
		 } else
		 { 
		  while(true)
		  {
		   x=PMinusMax-(PMinusMax-PMinusMin)*G4UniformRand();
		   y=G4UniformRand();
		   if( y < G4Pow::GetInstance()->powA(1.0-x/PMinusMax,3.) ) break;
		  }
		  PMinusNew = x;
		 };		 

		 TMinusNew = SqrtS - PMinusNew;

	         Qminus = Ptarget.minus() - TMinusNew;
		 TPlusNew = TargMassT2 / TMinusNew;
	         Qplus = Ptarget.plus() - TPlusNew;

		} else 
//-------------------------------------------------------------------------------
		{// The target will fragment, the projectile will saved.
		 TPlusMin = std::sqrt( TargMassT2 + PZcms2 ) - PZcms;
		 TPlusMax = SqrtS - ProjMassT;
		 sqrtTPlusMin = std::sqrt(TPlusMin); sqrtTPlusMax = std::sqrt(TPlusMax);

		 if( absPDGcode > 1000 ) 
		 { 
		  while(true)
		  {
		   x=TPlusMax-(TPlusMax-TPlusMin)*G4UniformRand();
		   y=G4UniformRand();
		   if( y < G4Pow::GetInstance()->powA(1.0-x/TPlusMax,3.) ) break;
		  }
		  TPlusNew = x;
		 } else if( (absPDGcode == 211) || (absPDGcode == 111) ) 
		 {
		  while(true)
		  {
		   x=sqrtTPlusMax-(sqrtTPlusMax-sqrtTPlusMin)*G4UniformRand();
		   y=G4UniformRand();
		   if( y < 1.0-0.7 * x/sqrtTPlusMax ) break;                  // 0.7 for Pi Vova
		  }
		  TPlusNew = sqr(x);
		 } else if( (absPDGcode == 321) || (absPDGcode == 311)  || 
			    (   PDGcode == 130) || (   PDGcode == 310) ) 
		 { // For K-mesons it must be found !!! Vova
		  while(true)
		  {
		   x=sqrtTPlusMax-(sqrtTPlusMax-sqrtTPlusMin)*G4UniformRand();
		   y=G4UniformRand();
		   if( y < 1.0-0.7 * x/sqrtTPlusMax ) break;
		  }
		 } else
		 { 
		  while(true)
		  {
		   x=TPlusMax-(TPlusMax-TPlusMin)*G4UniformRand();
		   y=G4UniformRand();
		   if( y < G4Pow::GetInstance()->powA(1.0-x/TPlusMax,3.) ) break;
		  }
		  TPlusNew = sqr(x);
		 };

		 PPlusNew = SqrtS - TPlusNew;

		 Qplus = PPlusNew - Pprojectile.plus();
		 PMinusNew = ProjMassT2 / PPlusNew;
		 Qminus = PMinusNew - Pprojectile.minus();
		}


	        Qmomentum.setPz( (Qplus - Qminus)/2 );
	        Qmomentum.setE(  (Qplus + Qminus)/2 );

  #ifdef debugQuarkExchange
     G4cout<<"ProjectileDiffraction (Pprojectile + Qmomentum).mag2()  Mprojectile2"<<G4endl;
     G4cout<<ProjectileDiffraction<<" "<<( Pprojectile + Qmomentum ).mag2()<<" "<< Mprojectile2<<G4endl;
     G4cout<<"TargetDiffraction     (Ptarget     - Qmomentum).mag2()  Mtarget2"<<G4endl;
     G4cout<<!ProjectileDiffraction<<" "<<( Ptarget    - Qmomentum ).mag2()<<" "<< Mtarget2<<G4endl;
  #endif

	} while ( ( ProjectileDiffraction&&( Pprojectile + Qmomentum ).mag2() <  Mprojectile2 ) ||
                  (!ProjectileDiffraction&&( Ptarget     - Qmomentum ).mag2() <  Mtarget2       )   ); 
                // Repeat the sampling because there was not any excitation


	Pprojectile += Qmomentum;

	Ptarget     -= Qmomentum;

	// Transform back and update SplitableHadron Participant.
	Pprojectile.transform(toLab);
	Ptarget.transform(toLab);

  #ifdef debugQuarkExchange
     G4cout << "Pprojectile  in Lab. " << Pprojectile << G4endl;
     G4cout << "Ptarget      in Lab. " << Ptarget     << G4endl;
     G4cout << "G4QuarkExchange: Projectile mass  " <<  Pprojectile.mag() << G4endl;
     G4cout << "G4QuarkExchange: Target mass      " <<  Ptarget.mag() << G4endl;
  #endif

	target->Set4Momentum(Ptarget);
	projectile->Set4Momentum(Pprojectile);

//=================================== Quark exchange ================================
	projectile->SplitUp();
	target->SplitUp();

	G4Parton* PrQuark = projectile->GetNextParton(); 
	G4Parton* TrQuark = target->GetNextParton(); 

	G4ParticleDefinition * Tmp = PrQuark->GetDefinition();
	PrQuark->SetDefinition(TrQuark->GetDefinition());
	TrQuark->SetDefinition(Tmp);

//============================================
	return true;
}


// --------- private methods ----------------------

G4ThreeVector G4QuarkExchange::GaussianPt(G4double widthSquare, G4double maxPtSquare) const
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





