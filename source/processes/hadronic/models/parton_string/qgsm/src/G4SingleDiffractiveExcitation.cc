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
// $Id: G4SingleDiffractiveExcitation.cc 106980 2017-10-31 09:02:49Z gcosmo $
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


G4SingleDiffractiveExcitation::G4SingleDiffractiveExcitation(G4double sigmaPt, G4double minextraMass,G4double x0mass)
: widthOfPtSquare(-2*sqr(sigmaPt)) , minExtraMass(minextraMass), minmass(x0mass)
{}

G4SingleDiffractiveExcitation::~G4SingleDiffractiveExcitation()
{}

G4bool G4SingleDiffractiveExcitation::
ExciteParticipants(G4VSplitableHadron *projectile, G4VSplitableHadron *target) const
{
  G4LorentzVector Pprojectile=projectile->Get4Momentum();
  G4double Mprojectile2=sqr(projectile->GetDefinition()->GetPDGMass() + minExtraMass);

  G4LorentzVector Ptarget=target->Get4Momentum();
  G4double Mtarget2=sqr(target->GetDefinition()->GetPDGMass() + minExtraMass);
  //G4cout << "E proj, target :" << Pprojectile.e() << ", " << Ptarget.e() << G4endl;

  G4bool KeepProjectile= G4UniformRand() > 0.5;

  // reset the min.mass of the non diffractive particle to its value, ( minus a bit for rounding...)
  if ( KeepProjectile )
  {
    //cout << " Projectile fix" << G4endl;
    Mprojectile2 = sqr(projectile->GetDefinition()->GetPDGMass() * (1-perCent) );
  } else {
    //cout << " Target fix" << G4endl;
    Mtarget2=sqr(target->GetDefinition()->GetPDGMass() * (1-perCent) );
  }

  // Transform momenta to cms and then rotate parallel to z axis;

  G4LorentzVector Psum;
  Psum=Pprojectile+Ptarget;

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

  //G4cout << "Pprojectile  be4 boost " << Pprojectile << G4endl;
  //G4cout << "Ptarget be4 boost : " << Ptarget << G4endl;

  G4LorentzRotation toLab(toCms.inverse());

  Pprojectile.transform(toCms);
  Ptarget.transform(toCms);

  G4LorentzVector Qmomentum;
  G4int whilecount=0;
  do {
    //  Generate pt

    G4double maxPtSquare=sqr(Ptarget.pz());
    if (whilecount++ >= 500 && (whilecount%100)==0)
      //G4cout << "G4SingleDiffractiveExcitation::ExciteParticipants possibly looping"
      //       << ", loop count/ maxPtSquare : "
      //       << whilecount << " / " << maxPtSquare << G4endl;

    if (whilecount > 1000 )
    {
      Qmomentum=G4LorentzVector(0.,0.,0.,0.);
      //G4cout << "G4SingleDiffractiveExcitation::ExciteParticipants: Aborting loop!" << G4endl;
      return false; 	  //  Ignore this interaction
    }
    Qmomentum=G4LorentzVector(GaussianPt(widthOfPtSquare,maxPtSquare),0);

    //  Momentum transfer
    G4double Xmin = minmass / ( Pprojectile.e() + Ptarget.e() );
    G4double Xmax=1.;
    G4double Xplus =ChooseX(Xmin,Xmax);
    G4double Xminus=ChooseX(Xmin,Xmax);

    G4double pt2=G4ThreeVector(Qmomentum.vect()).mag2();
    G4double Qplus =-1 * pt2 / Xminus/Ptarget.minus();
    G4double Qminus=     pt2 / Xplus /Pprojectile.plus();

    if ( KeepProjectile )
    {
      Qminus = (sqr(projectile->GetDefinition()->GetPDGMass()) + pt2 )
	       / (Pprojectile.plus() + Qplus ) -  Pprojectile.minus();
    } else {
      Qplus = Ptarget.plus() - (sqr(target->GetDefinition()->GetPDGMass()) + pt2 )
		               / (Ptarget.minus() - Qminus );
    }			

    Qmomentum.setPz( (Qplus-Qminus)/2 );
    Qmomentum.setE(  (Qplus+Qminus)/2 );

    //G4cout << "Qplus / Qminus " << Qplus << " / " << Qminus<<G4endl;
    //G4cout << "pt2 " << pt2 << G4endl;
    //G4cout << "Qmomentum " << Qmomentum << G4endl;
    //G4cout << " Masses (P/T) : " << (Pprojectile+Qmomentum).mag() <<
    //	        " / " << (Ptarget-Qmomentum).mag() << G4endl;

  } while (  (Ptarget-Qmomentum).mag2() <= Mtarget2  /* Loop checking, 26.10.2015, A.Ribon */
	     || (Pprojectile+Qmomentum).mag2() <= Mprojectile2
	     || (Ptarget-Qmomentum).e() < 0.
	     || (Pprojectile+Qmomentum).e() < 0. );

  //G4double Ecms=Pprojectile.e() + Ptarget.e();

  Pprojectile += Qmomentum;
  Ptarget     -= Qmomentum;

  //G4cout << "Pprojectile.e()  : " << Pprojectile.e() << G4endl;
  //G4cout << "Ptarget.e()      : " << Ptarget.e() << G4endl;
  //G4cout << "end event_______________________________________________"<<G4endl;
  //G4cout << "Pprojectile with Q : " << Pprojectile << G4endl;
  //G4cout << "Ptarget with Q : " << Ptarget << G4endl;
  //G4cout << "Projectile back: " << toLab * Pprojectile << G4endl;
  //G4cout << "Target back: " << toLab * Ptarget << G4endl;

  // Transform back and update SplitableHadron Participant.
  Pprojectile.transform(toLab);
  Ptarget.transform(toLab);

  //G4cout << "G4SingleDiffractiveExcitation- Target mass      " <<  Ptarget.mag() << G4endl;
  //G4cout << "G4SingleDiffractiveExcitation- Projectile mass  " <<  Pprojectile.mag() << G4endl;

  target->Set4Momentum(Ptarget);
  projectile->Set4Momentum(Pprojectile);

  return true;
}


// --------- private methods ----------------------

G4double G4SingleDiffractiveExcitation::ChooseX(G4double Xmin, G4double Xmax) const
{
  // choose an x between Xmin and Xmax with P(x) ~ 1/x
  //  to be improved...

  G4double range=Xmax-Xmin;

  if ( Xmin <= 0. || range <=0. ) 
  {
    G4cout << " Xmin, range : " << Xmin << " , " << range << G4endl;
    throw G4HadronicException(__FILE__, __LINE__, "G4SingleDiffractiveExcitation::ChooseX : Invalid arguments ");
  }

  G4double x;
  do {
    x=Xmin + G4UniformRand() * range;
  } while ( Xmin/x < G4UniformRand() );  /* Loop checking, 26.10.2015, A.Ribon */

  //cout << "DiffractiveX "<<x<<G4endl;
  return x;
}

G4ThreeVector G4SingleDiffractiveExcitation::GaussianPt(G4double widthSquare, G4double maxPtSquare) const
{
  //  @@ this method is used in FTFModel as well. Should go somewhere common!

  G4double pt2;

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = -1;
  do {
    pt2=widthSquare * G4Log( G4UniformRand() );
  } while ( ( pt2 > maxPtSquare) && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 26.10.2015, A.Ribon */
  if ( loopCounter >= maxNumberOfLoops ) pt2 = 0.0;

  pt2=std::sqrt(pt2);

  G4double phi=G4UniformRand() * twopi;

  return G4ThreeVector (pt2*std::cos(phi), pt2*std::sin(phi), 0.);    
}

