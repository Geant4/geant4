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
	
	
