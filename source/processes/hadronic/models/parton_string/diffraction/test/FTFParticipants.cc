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
// $Id: FTFParticipants.cc,v 1.1 2003-10-08 13:48:52 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 file
//
//
//             by Gunter Folger, June 1998.
//       class exercising G4FTFParticipants class.
// ------------------------------------------------------------

#include "G4FTFParticipants.hh"

#include "G4StableIsotopes.hh"
#include "G4ReactionProduct.hh"
#include "G4PionPlus.hh"

// Solve templates.
// Code in g4templates.hh controlled by the macro G4_SOLVE_TEMPLATES
#include "g4templates.hh"
#ifdef G4_SOLVE_TEMPLATES
#include "G4InteractionContent.hh"
template class G4RWTPtrOrderedVector<G4InteractionContent>;
#include "G4VSplitableHadron.hh"
template class G4RWTPtrOrderedVector<G4VSplitableHadron>;
#endif

int main()
{

	G4FTFParticipants participant;
	
	G4ReactionProduct projectile(G4PionPlus::PionPlus());
	G4cout << " pion - mass : " << G4PionPlus::PionPlus()->GetPDGMass() /GeV
	<< G4endl;
	G4cout  << " proj - mass : " << projectile.GetMass() / GeV << G4endl;
	G4ThreeVector momentum(100*GeV,100*GeV,100*GeV);
	projectile.SetMomentum(momentum);
	projectile.SetTotalEnergy(sqrt(sqr(projectile.GetMass())+momentum.mag2()));
	
	G4cout  << " proj - Etot : " << projectile.GetTotalEnergy() / GeV << G4endl;

	G4StableIsotopes theIso;	

//	for (int Z=1; Z<93; Z++ )
	for (int Z=18; Z<19; Z++ )
	{
	   G4cout <<G4endl<< G4endl<< "new Element " ;
	   G4cout << theIso.GetName(Z) << G4endl;
	   G4cout << "         with " << theIso.GetNumberOfIsotopes(Z) ;
	   G4cout << " Isotopes"<< G4endl;
	   for (G4int iso=0; iso<theIso.GetNumberOfIsotopes(Z); iso++)
	   {
	       G4double massnumber=
	          theIso.GetIsotopeNucleonCount(theIso.GetFirstIsotope(Z)+iso);
	       participant.Init(massnumber,(G4double) Z);
	       G4cout << "Charge" << Z << G4endl;
	       G4cout << "Mass number" << massnumber << G4endl;
	      
	      participant.BuildInteractions(projectile);
	   }   
	}
	return 0;
	
} 
	
	
