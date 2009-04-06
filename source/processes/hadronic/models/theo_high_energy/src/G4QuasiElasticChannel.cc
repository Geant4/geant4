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
// $Id: G4QuasiElasticChannel.cc,v 1.6 2009-04-06 09:35:46 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Author : Gunter Folger March 2007
// Modified by Mikhail Kossov. Apr2009, E/M conservation: ResidualNucleus is added (ResNuc)
// Class Description
// Final state production model for theoretical models of hadron inelastic
// quasi elastic scattering in geant4;
// Class Description - End
//
// Modified:

#include "G4QuasiElasticChannel.hh"

#include "G4Fancy3DNucleus.hh"


#include "G4HadTmpUtil.hh"		//lrint

//#define debug_scatter

G4QuasiElasticChannel::G4QuasiElasticChannel()
{
	theQuasiElastic=G4QuasiFreeRatios::GetPointer();
}

G4QuasiElasticChannel::~G4QuasiElasticChannel()
{}

G4double G4QuasiElasticChannel::GetFraction(G4Nucleus &theNucleus,
				const G4DynamicParticle & thePrimary)
{
	std::pair<G4double,G4double> ratios;
	ratios=theQuasiElastic->GetRatios(thePrimary.GetTotalMomentum(),
			thePrimary.GetDefinition()->GetPDGEncoding(),
			G4lrint(theNucleus.GetZ()),
			G4lrint(theNucleus.GetN()-theNucleus.GetZ()));
#ifdef debug_getFraction			
	G4cout << "G4QuasiElasticChannel::ratios " << ratios.first << " x " <<ratios.second
	       << "  = " << ratios.first*ratios.second << G4endl;
#endif
	       
	//return ratios.first*ratios.second;
	//return 0.; // Switch off quasielastic (temporary) M.K.
	return 1.; // Only quasielastic (temporary) M.K.
}

G4KineticTrackVector * G4QuasiElasticChannel::Scatter(G4Nucleus &theNucleus,
					const G4DynamicParticle & thePrimary)
{
	
	
	G4int A=G4lrint(theNucleus.GetN());
	G4int Z=G4lrint(theNucleus.GetZ());                                 // M.K. ResNuc
//   build Nucleus and choose random nucleon to scatter with
	the3DNucleus.Init(theNucleus.GetN(),theNucleus.GetZ());
	const std::vector<G4Nucleon *> nucleons=the3DNucleus.GetNucleons();
      G4double targetNucleusMass=the3DNucleus.GetMass();                  // M.K. ResNuc
	G4LorentzVector targetNucleus4Mom(0.,0.,0.,targetNucleusMass);      // M.K. ResNuc
      G4int index;
      do {
	   index=G4lrint((A-1)*G4UniformRand());
	   
	} while (index < 0 || index >= static_cast<G4int>(nucleons.size()));
	    
//	G4double an=G4UniformRand()*A;
	G4ParticleDefinition * pDef= nucleons[index]->GetDefinition();
	G4int resA=A-1;                                                     // M.K. ResNuc
      G4int resZ=Z-static_cast<int>(pDef->GetPDGCharge());                // M.K. ResNuc
      G4ParticleDefinition* resDef=G4ParticleTable::GetParticleTable()
                                             ->FindIon(resZ,resA,0,resZ); // M.K. ResNuc
      G4double residualNucleusMass=resDef->GetPDGMass();                  // M.K. ResNuc
#ifdef debug_scatter
	G4cout << " neutron - proton? A, Z, an, pdg" <<" "
	       << A <<" "<< Z << " "<<an <<" " << pDef->GetParticleName()<< G4endl;
#endif
//	G4LorentzVector pNucleon(G4ThreeVector(0,0,0),pDef->GetPDGMass());
	G4LorentzVector pNucleon=nucleons[index]->Get4Momentum();
      G4double residualNucleusEnergy=std::sqrt(residualNucleusMass*residualNucleusMass+
                                               pNucleon.vect().mag2());   // M.K. ResNuc
      pNucleon.setE(targetNucleusMass-residualNucleusEnergy);             // M.K. ResNuc
      G4LorentzVector residualNucleus4Mom=targetNucleus4Mom-pNucleon;     // M.K. ResNuc
	
	std::pair<G4LorentzVector,G4LorentzVector> result;
	
	result=theQuasiElastic->Scatter(pDef->GetPDGEncoding(),pNucleon,
	                       thePrimary.GetDefinition()->GetPDGEncoding(),
			        thePrimary.Get4Momentum());
      G4LorentzVector scatteredHadron4Mom=result.second;                  // M.K. ResNuc
	if (result.first.e() <= 0.)
	{
        //G4cout << "Warning - G4QuasiElasticChannel::Scatter no scattering" << G4endl;
	  //return 0;       //no scatter
        G4LorentzVector scatteredHadron4Mom=thePrimary.Get4Momentum();   // M.K. ResNuc
        residualNucleus4Mom=G4LorentzVector(0.,0.,0.,targetNucleusMass); // M.K. ResNuc
        resDef=G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);    // M.K. ResNuc
	}

#ifdef debug_scatter
        G4LorentzVector EpConservation=pNucleon+thePrimary.Get4Momentum() 
	                                - result.first - result.second;
	if (   (EpConservation.vect().mag2() > 0.01*MeV*MeV )
	    || 	(std::abs(EpConservation.e()) > 0.1 * MeV ) ) 
	{
	    G4cout << "Warning - G4QuasiElasticChannel::Scatter E-p non conservation : "
	            << EpConservation << G4endl;
	}    
#endif

	G4KineticTrackVector * ktv;
	ktv=new G4KineticTrackVector();
	G4KineticTrack * sPrim=new G4KineticTrack(thePrimary.GetDefinition(),
                                                0.,G4ThreeVector(0), scatteredHadron4Mom);
      ktv->push_back(sPrim);
	if (result.first.e() > 0.)
	{
        G4KineticTrack * sNuc=new G4KineticTrack(pDef, 0.,G4ThreeVector(0), result.first);
        ktv->push_back(sNuc);
      }
	G4KineticTrack * rNuc=new G4KineticTrack(resDef,
                               0.,G4ThreeVector(0), residualNucleus4Mom); // M.K. ResNuc
	ktv->push_back(rNuc);                                               // M.K. ResNuc
	
#ifdef debug_scatter
	G4cout << " scattered Nucleon : " << result.first << " mass " <<result.first.mag() << G4endl;
	G4cout << " scattered Project : " << result.second << " mass " <<result.second.mag() << G4endl;
#endif
	return ktv;			       
}
