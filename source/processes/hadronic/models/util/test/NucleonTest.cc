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
// $Id: NucleonTest.cc,v 1.1 2003-10-08 13:48:24 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 file
//
//
//             by Gunter Folger, June 1998.
//       class exercising G4Nucleon class.
// ------------------------------------------------------------

#include "G4Nucleon.hh"

// Solve templates.
// Code in g4templates.hh controlled by the macro G4_SOLVE_TEMPLATES
#include "g4templates.hh"


int main()
{

	G4Nucleon nucleon;
	G4ThreeVector vector;
        G4ThreeVector beta(0.51,-0.35,0.17); 
        G4LorentzVector init(0,0,400.,1000.);
	
//                             mass**2   + P**2
	G4double energy=init.e();   //sqrt( init.mag2() +beta.dot(beta)*init.vect().mag2() );
        G4cout << " Energy : "  << energy << G4endl;
	G4LorentzVector
	Lboost(-1/sqrt(1-beta.mag2())*beta,sqrt(1 + beta.mag2()/(1-beta.mag2())));
	vector= G4ThreeVector(1.,2.,5.);
	
	nucleon.SetPosition(vector);
	
	G4ThreeVector pos = nucleon.GetPosition();
	G4cout << " pos.x :" << pos.x() << G4endl;
	G4cout << " pos.y :" << pos.y() << G4endl;
	G4cout << " pos.z :" << pos.z() << G4endl;
	
	
	nucleon.SetMomentum(Lboost);

	G4cout << "Momentum for boost: " <<
	      nucleon.GetMomentum().px() << " , " <<
	      nucleon.GetMomentum().py() << " , " <<
	      nucleon.GetMomentum().pz() << G4endl <<
	      "Energy/mass : " << nucleon.GetMomentum().e() <<  " , " <<
	      nucleon.GetMomentum().mag() << G4endl << G4endl;
	
	
	nucleon.SetMomentum(init);

	G4cout << "Momentum before boost: " <<
	      nucleon.GetMomentum().px() << " , " <<
	      nucleon.GetMomentum().py() << " , " <<
	      nucleon.GetMomentum().pz() << G4endl <<
	      "Energy/mass : " << nucleon.GetMomentum().e() <<  " , " <<
	      nucleon.GetMomentum().mag() << G4endl << G4endl;
	
	nucleon.Boost(beta);
	
	G4cout << "Momentum after boost: " <<
	      nucleon.GetMomentum().px() << " , " <<
	      nucleon.GetMomentum().py() << " , " <<
	      nucleon.GetMomentum().pz() << G4endl <<
	      "Energy/mass : " << nucleon.GetMomentum().e() <<  " , " <<
	      nucleon.GetMomentum().mag() << G4endl << G4endl;
	      
	nucleon.SetMomentum(init);
	
	nucleon.Boost(Lboost);
	
	G4cout << "Momentum after Lboost: " <<
	      nucleon.GetMomentum().px() << " , " <<
	      nucleon.GetMomentum().py() << " , " <<
	      nucleon.GetMomentum().pz() << G4endl <<
	      "Energy/mass : " << nucleon.GetMomentum().e() << " , " <<
	      nucleon.GetMomentum().mag() << G4endl << G4endl;
	      
	      
	      
	      
} 
	
	
