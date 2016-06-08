// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DiffractiveExcitation.cc,v 1.7 1999/12/15 17:51:26 gcosmo Exp $
// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4DiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// ------------------------------------------------------------


#include "globals.hh"
#include "Randomize.hh"

#include "G4DiffractiveExcitation.hh"
#include "G4LorentzRotation.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"
//#include "G4ios.hh"

G4DiffractiveExcitation::G4DiffractiveExcitation(G4double sigmaPt, G4double minextraMass,G4double x0mass)
:
widthOfPtSquare(-2*sqr(sigmaPt)) , minExtraMass(minextraMass),
minmass(x0mass)
{
}

G4bool G4DiffractiveExcitation::
  ExciteParticipants(G4VSplitableHadron *projectile, G4VSplitableHadron *target) const
{

	   G4LorentzVector Pprojectile=projectile->Get4Momentum();
	   G4double Mprojectile2=sqr(projectile->GetDefinition()->GetPDGMass() + minExtraMass);

  	   G4LorentzVector Ptarget=target->Get4Momentum();
   	   G4double Mtarget2=sqr(target->GetDefinition()->GetPDGMass() + minExtraMass);
//	     G4cout << "E proj, target :" << Pprojectile.e() << ", " <<
//					    Ptarget.e() << G4endl;

// Transform momenta to cms and then rotate parallel to z axis;

	   G4LorentzVector Psum;
	   Psum=Pprojectile+Ptarget;

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

//	   G4cout << "Pprojectile  be4 boost " << Pprojectile << G4endl;
//	   G4cout << "Ptarget be4 boost : " << Ptarget << G4endl;
	


	   G4LorentzRotation toLab(toCms.inverse());

	   Pprojectile.transform(toCms);
	   Ptarget.transform(toCms);

//	   G4cout << "Pprojectile aft boost : " << Pprojectile << G4endl;
//	   G4cout << "Ptarget aft boost : " << Ptarget << G4endl;
//	   G4cout << "cms aft boost : " << (Pprojectile+ Ptarget) << G4endl;

//	   G4cout << " Projectile Xplus / Xminus : " << 
//	   	Pprojectile.plus() << " / " << Pprojectile.minus() << G4endl;
//	   G4cout << " Target Xplus / Xminus : " << 
//	   	Ptarget.plus() << " / " << Ptarget.minus() << G4endl;

	   G4LorentzVector Qmomentum;
	   G4int whilecount=0;
	   do {
//  Generate pt		

	       G4double maxPtSquare=sqr(Ptarget.pz());
	       if (whilecount++ >= 500 && (whilecount%100)==0) 
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : " 
//           	 << whilecount << " / " << maxPtSquare << G4endl;
               if (whilecount > 1000 ) 
               {
           	   Qmomentum=0;
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants: Aborting loop!" << G4endl;
	   	 return false; 	  //  Ignore this interaction 
               }
	       Qmomentum=G4LorentzVector(GaussianPt(widthOfPtSquare,maxPtSquare),0);

//	       G4cout << "generated Pt " << Qmomentum << G4endl;
//	       G4cout << "Pprojectile with pt : " << Pprojectile+Qmomentum << G4endl;
//	       G4cout << "Ptarget with pt : " << Ptarget-Qmomentum << G4endl;

//  Momentum transfer
	       G4double Xmin = minmass / ( Pprojectile.e() + Ptarget.e() );
	       G4double Xmax=1.;
	       G4double Xplus =ChooseX(Xmin,Xmax);
	       G4double Xminus=ChooseX(Xmin,Xmax);

//	       G4cout << " X-plus  " << Xplus << G4endl;
//	       G4cout << " X-minus " << Xminus << G4endl;
	       
	       G4double pt2=G4ThreeVector(Qmomentum).mag2();
	       G4double Qplus =-1 * pt2 / Xminus/Ptarget.minus();
	       G4double Qminus=     pt2 / Xplus /Pprojectile.plus();

	       Qmomentum.setPz( (Qplus-Qminus)/2 );
	       Qmomentum.setE(  (Qplus+Qminus)/2 );

//	     G4cout << "Qplus / Qminus " << Qplus << " / " << Qminus<<G4endl;
//	     G4cout << "pt2" << pt2 << G4endl;
//	     G4cout << "Qmomentum " << Qmomentum << G4endl;
//	     G4cout << " Masses (P/T) : " << (Pprojectile+Qmomentum).mag() <<
//	   		       " / " << (Ptarget-Qmomentum).mag() << G4endl;

	   } while ( (Pprojectile+Qmomentum).mag2() <= Mprojectile2 ||
	   	     (Ptarget-Qmomentum).mag2()     <= Mtarget2 );
	   			   
	   Pprojectile += Qmomentum;
	   Ptarget     -= Qmomentum;

//	   G4cout << "Pprojectile with Q : " << Pprojectile << G4endl;
//	   G4cout << "Ptarget with Q : " << Ptarget << G4endl;
	
	
	
//	   G4cout << "Projectile back: " << toLab * Pprojectile << G4endl;
//	   G4cout << "Target back: " << toLab * Ptarget << G4endl;

// Transform back and update SplitableHadron Participant.
	   Pprojectile.transform(toLab);
	   Ptarget.transform(toLab);

//	   G4cout << "Target	 mass  " <<  Ptarget.mag() << G4endl;
	   		   
	   target->Set4Momentum(Ptarget);
// 	
	   

//	   G4cout << "Projectile mass  " <<  Pprojectile.mag() << G4endl;

	   projectile->Set4Momentum(Pprojectile);
	
	
	return true;
}


G4ExcitedString * G4DiffractiveExcitation::
           String(G4VSplitableHadron * hadron, G4bool isProjectile) const
{
	hadron->SplitUp();
	G4Parton *start= hadron->GetNextParton();
	if ( start==NULL) 
	{ G4cout << " G4FTFModel::String() Error:No start parton found"<< G4endl;
	  return NULL;
	}
	G4Parton *end  = hadron->GetNextParton();
	if ( end==NULL) 
	{ G4cout << " G4FTFModel::String() Error:No end parton found"<< G4endl;
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
	G4double ptSquared= hadron->Get4Momentum().perp2();
	G4double transverseMassSquared= hadron->Get4Momentum().plus()
				    *	hadron->Get4Momentum().minus();


	G4double maxAvailMomentumSquared=
		 sqr( sqrt(transverseMassSquared) - sqrt(ptSquared) );

	G4ThreeVector pt=GaussianPt(widthOfPtSquare,maxAvailMomentumSquared);

	G4LorentzVector Pstart(G4LorentzVector(pt,0.));
	G4LorentzVector Pend;
	Pend.setPx(hadron->Get4Momentum().px() - pt.x());
	Pend.setPy(hadron->Get4Momentum().py() - pt.y());

	G4double tm1=hadron->Get4Momentum().minus() +
	  ( Pend.perp2()-Pstart.perp2() ) / hadron->Get4Momentum().plus();

	G4double tm2= sqrt( G4std::max(0., sqr(tm1) -
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
	
#ifdef G4_FTFDEBUG
	G4cout << " generated string flavors          " << start->GetPDGcode() << " / " << end->GetPDGcode() << G4endl;
	G4cout << " generated string momenta:   quark " << start->Get4Momentum() << "mass : " <<start->Get4Momentum().mag()<< G4endl;
	G4cout << " generated string momenta: Diquark " << end ->Get4Momentum() << "mass : " <<end->Get4Momentum().mag()<< G4endl;
	G4cout << " sum of ends                       " << Pstart+Pend << G4endl;
	G4cout << " Original                          " << hadron->Get4Momentum() << G4endl;
#endif
	
	return string;	
}


// --------- private methods ----------------------

G4double G4DiffractiveExcitation::ChooseX(G4double Xmin, G4double Xmax) const
{
// choose an x between Xmin and Xmax with P(x) ~ 1/x

//  to be improved...

	G4double range=Xmax-Xmin;
	
	if ( Xmin <= 0. || range <=0. ) 
	{
		G4cout << " Xmin, range : " << Xmin << " , " << range << G4endl;
		G4Exception("G4DiffractiveExcitation::ChooseX : Invalid arguments ");
	}

	G4double x;
	do {
	    x=Xmin + G4UniformRand() * range;
	}  while ( Xmin/x < G4UniformRand() );

//debug-hpw	cout << "DiffractiveX "<<x<<G4endl;
	return x;
}

G4ThreeVector G4DiffractiveExcitation::GaussianPt(G4double widthSquare, G4double maxPtSquare) const
{            //  @@ this method is used in FTFModel as well. Should go somewhere common!
	
	G4double pt2;

	do {
	    pt2=widthSquare * log( G4UniformRand() );
	} while ( pt2 > maxPtSquare);
	
	pt2=sqrt(pt2);
	
	G4double phi=G4UniformRand() * twopi;
	
	return G4ThreeVector (pt2*cos(phi), pt2*sin(phi), 0.);    
}

G4DiffractiveExcitation::G4DiffractiveExcitation(const G4DiffractiveExcitation &right)
:
widthOfPtSquare(0) , minExtraMass(0),
minmass(0)
{
	G4Exception("G4DiffractiveExcitation copy contructor not meant to be called");
}


G4DiffractiveExcitation::~G4DiffractiveExcitation()
{
}


const G4DiffractiveExcitation & G4DiffractiveExcitation::operator=(const G4DiffractiveExcitation &right)
{
	G4Exception("G4DiffractiveExcitation = operator meant to be called");
	return *this;
}


int G4DiffractiveExcitation::operator==(const G4DiffractiveExcitation &right) const
{
	G4Exception("G4DiffractiveExcitation == operator meant to be called");
	return false;
}

int G4DiffractiveExcitation::operator!=(const G4DiffractiveExcitation &right) const
{
	G4Exception("G4DiffractiveExcitation != operator meant to be called");
	return true;
}





