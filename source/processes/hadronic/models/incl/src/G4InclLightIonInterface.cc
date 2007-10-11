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
// $Id: G4InclLightIonInterface.cc,v 1.3 2007-10-11 08:15:37 gcosmo Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "G4InclLightIonInterface.hh"
#include "math.h"

G4InclLightIonInterface::G4InclLightIonInterface()
{
  hazard = new G4Hazard();
  hazard->ial = 37723;
  hazard->igraine[0] = 3997;
  hazard->igraine[1] = 15573;
  hazard->igraine[2] = 9971;
  hazard->igraine[3] = 9821; 
  hazard->igraine[4] = 99233; 
  hazard->igraine[5] = 11167; 
  hazard->igraine[6] = 12399;
  hazard->igraine[7] = 11321; 
  hazard->igraine[8] = 9825;
  hazard->igraine[9] = 2587; 
  hazard->igraine[10] = 1775;
  hazard->igraine[11] = 56799; 
  hazard->igraine[12] = 1156;
  //  hazard->igraine[13] = 11207;
  hazard->igraine[13] = 38957; 
  hazard->igraine[14] = 35779; 
  hazard->igraine[15] = 10055; 
  hazard->igraine[16] = 76533; 
  hazard->igraine[17] = 33759;
  hazard->igraine[18] = 13227;

}

G4InclLightIonInterface::~G4InclLightIonInterface()
{

}

// For singleton implementation:

G4HadFinalState* G4InclLightIonInterface::ApplyYourself(const G4HadProjectile& aTrack, 
			       G4Nucleus& theNucleus)
{
  G4VarNtp *varntp = new G4VarNtp();
  G4Calincl *calincl = new G4Calincl();
  G4Ws *ws = new G4Ws();
  G4Mat *mat = new G4Mat();

  G4Incl *incl = new G4Incl(hazard, calincl, ws, mat, varntp);

  G4int maxTries = 200;
  G4int tries = 0;

  G4int particleI, n = 0;

  G4int bulletType = 0;

  // Print diagnostic messages. 0 = silent, 1 and 2 = verbose
  //  verboseLevel = 3;

  // Increase the event number:
  eventNumber++;

  if (verboseLevel > 0) {
    G4cout << " >>> G4InclLightIonInterface::ApplyYourself called" << G4endl;
  }

  if(verboseLevel > 1) {
    G4cout <<"G4InclLightIonInterface: Now processing INCL4 event number:" << eventNumber << G4endl;
  }

  // INCL4 needs the energy in units MeV
  G4double bulletE = aTrack.GetKineticEnergy() / MeV;

  G4double targetA = theNucleus.GetN();
  G4double targetZ = theNucleus.GetZ();

  G4double eKin;
  G4double momx, momy, momz;
  G4DynamicParticle *cascadeParticle = 0;
  G4ParticleDefinition *aParticleDefinition = 0;

  // INCL assumes the projectile particle is going in the direction of
  // the Z-axis. Here we construct proper rotation to convert the
  // momentum vectors of the outcoming particles to the original
  // coordinate system.
  G4LorentzVector projectileMomentum = aTrack.Get4Momentum();
  G4LorentzRotation toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4LorentzRotation toLabFrame = toZ.inverse();

  // Make sure the output data structure is clean.
  theResult.Clear();

  // Map Geant4 particle types to corresponding INCL4 types.
  enum bulletParticleType {nucleus = 0, proton = 1, neutron = 2, pionPlus = 3, pionZero = 4, 
                           pionMinus = 5, deuteron = 6, triton = 7, he3 = 8, he4 = 9};

  // Coding particles for use with INCL4 and ABLA 
  if (aTrack.GetDefinition() == G4Deuteron::Deuteron()   ) bulletType = deuteron;
  if (aTrack.GetDefinition() == G4Triton::Triton()       ) bulletType = triton;
  if (aTrack.GetDefinition() == G4He3::He3()             ) bulletType = he3;
  if (aTrack.GetDefinition() == G4Alpha::Alpha()         ) bulletType = he4;

  // Clean up FINPUT
  for(int i = 0; i < 15; i++) {
    calincl->f[i] = 0.0;
  }
  // End of initial clean up.

  // Check wheter the input is acceptable. This will contain more tests in the future. 
  if((bulletType != 0) && ((targetA != 1) && (targetZ != 1))) {
    // Set input parameters:

    // Target nucleus:
    // Mass number:
    calincl->f[0] = targetA;

    // Charge number:
    calincl->f[1] = targetZ;

    // Bullet:
    // Bullet type (1.00 = proton):
    calincl->f[6] = bulletType;
    // Bullet energy:
    calincl->f[2] = bulletE;

    // Run settings:
    // Time scaling:
    calincl->f[5] = 1.0;

    // Nuclear potential:
    calincl->f[4] = 45.0;

    // nosurf:
    ws->nosurf = -2;

    // XFOISA
    ws->xfoisa = 8;

    // NPAULSTR
    ws->npaulstr = 0;

    // Print event
    //    int  ievtest = 0; // deprecated in release alpha-2c

    // icoup variable is deprecated
    calincl->icoup = 1;
    
    tries = 0;

    varntp->ntrack = 0;
    G4int round = 0;

    // Loop until we produce real cascade
    mat->nbmat = 1;
    mat->amat[0] = int(calincl->f[0]);
    mat->zmat[0] = int(calincl->f[1]);

    // Initialize INCL
    incl->initIncl(true);

    while((varntp->ntrack <= 0) && (tries < maxTries)) {
      if(verboseLevel > 1) {
        round++;
        G4cout <<"Trying INCL4+ABLA. Round : " << round << G4endl; 
      }
      incl->processEventIncl();
      tries++;
      if(verboseLevel > 1) {
        G4cout <<"ntrack : " <<  varntp->ntrack <<G4endl;
      }
    }

    if(verboseLevel > 0) {
      // Diagnostic output
      G4cout <<"G4Incl4AblaCascadeInterface: Bullet type: " << bulletType << G4endl;
      G4cout <<"G4Incl4AblaCascadeInterface: Bullet energy: " << bulletE << " MeV" << G4endl;

      G4cout <<"G4Incl4AblaCascadeInterface: Target A:  " << targetA << G4endl;
      G4cout <<"G4Incl4AblaCascadeInterface: Target Z:  " << targetZ << G4endl;

      if(verboseLevel > 3) {
        diagdata <<"G4Incl4AblaCascadeInterface: Bullet type: " << bulletType << G4endl;
        diagdata <<"G4Incl4AblaCascadeInterface: Bullet energy: " << bulletE << " MeV" << G4endl;
        
        diagdata <<"G4Incl4AblaCascadeInterface: Target A:  " << targetA << G4endl;
        diagdata <<"G4Incl4AblaCascadeInterface: Target Z:  " << targetZ << G4endl;
      }

      for(particleI = 0; particleI < varntp->ntrack; particleI++) {
        G4cout << n << " " << calincl->f[6]  << " " << calincl->f[2] << " ";
        G4cout << varntp->massini << " " << varntp->mzini << " ";
        G4cout << varntp->exini << " " << varntp->mulncasc << " " << varntp->mulnevap << " " << varntp->mulntot << " ";
        G4cout << varntp->bimpact << " " << varntp->jremn << " " << varntp->kfis << " " << varntp->estfis << " ";
        G4cout << varntp->izfis << " " << varntp->iafis << " " << varntp->ntrack << " " << varntp->itypcasc[particleI] << " ";
        G4cout << varntp->avv[particleI] << " " << varntp->zvv[particleI] << " " << varntp->enerj[particleI] << " ";
        G4cout << varntp->plab[particleI] << " " << varntp->tetlab[particleI] << " " << varntp->philab[particleI] << G4endl;
        // For diagnostic output
        if(verboseLevel > 2) {
          diagdata << n << " " << calincl->f[6]  << " " << calincl->f[2] << " ";
          diagdata << varntp->massini << " " << varntp->mzini << " ";
          diagdata << varntp->exini << " " << varntp->mulncasc << " " << varntp->mulnevap << " " << varntp->mulntot << " ";
          diagdata << varntp->bimpact << " " << varntp->jremn << " " << varntp->kfis << " " << varntp->estfis << " ";
          diagdata << varntp->izfis << " " << varntp->iafis << " " << varntp->ntrack << " " << varntp->itypcasc[particleI] << " ";
          diagdata << varntp->avv[particleI] << " " << varntp->zvv[particleI] << " " << varntp->enerj[particleI] << " ";
          diagdata << varntp->plab[particleI] << " " << varntp->tetlab[particleI] << " " << varntp->philab[particleI] << G4endl;
        }
      }
    }

    // Check whether a valid cascade was produced. If not return the
    // original bullet particle with the same momentum.
    if(varntp->ntrack <= 0) {
      if(verboseLevel > 0) {
        G4cout <<"G4InclLightIonInterface: No cascade. Returning original particle with original momentum." << G4endl;
      }

      theResult.SetStatusChange(stopAndKill);
      
      if(bulletType == proton) {
        aParticleDefinition = G4Proton::ProtonDefinition();
      }
      if(bulletType == neutron) {
        aParticleDefinition = G4Neutron::NeutronDefinition();
      }
      if(bulletType == pionPlus) {
        aParticleDefinition = G4PionPlus::PionPlusDefinition();
      }
      if(bulletType == pionZero) {
        aParticleDefinition = G4PionZero::PionZeroDefinition();
      }
      if(bulletType == pionMinus) {
        aParticleDefinition = G4PionMinus::PionMinusDefinition();
      }
      if(bulletType == deuteron) {
	aParticleDefinition = G4Deuteron::DeuteronDefinition();
      }
      if(bulletType == triton) {
	aParticleDefinition = G4Triton::TritonDefinition();
      }
      if(bulletType == he3) {
	aParticleDefinition = G4He3::He3Definition();
      }
      if(bulletType == he4) {
	aParticleDefinition = G4Alpha::AlphaDefinition();
      }

      cascadeParticle = new G4DynamicParticle();
      cascadeParticle->SetDefinition(aParticleDefinition);
      cascadeParticle->Set4Momentum(aTrack.Get4Momentum());
      theResult.AddSecondary(cascadeParticle);
    }

    // Convert INCL4 output to Geant4 compatible data structures.
    // Elementary particles are converted to G4DynamicParticle.
    theResult.SetStatusChange(stopAndKill);
    // Loop through the INCL4+ABLA output.
    for(particleI = 0; particleI < varntp->ntrack; particleI++) {
      // Get energy/momentum and construct momentum vector:
      // In INCL4 coordinates!
      momx = varntp->plab[particleI]*std::cos(varntp->tetlab[particleI]*M_PI/180.0)*std::sin(varntp->philab[particleI]*M_PI/180.0)*MeV;
      momy = varntp->plab[particleI]*std::sin(varntp->tetlab[particleI]*M_PI/180.0)*std::sin(varntp->philab[particleI]*M_PI/180.0)*MeV;
      momz = varntp->plab[particleI]*std::cos(varntp->tetlab[particleI]*M_PI/180.0)*MeV;

      eKin = varntp->enerj[particleI] * MeV;

      //      eExcit = varntp->

      if(verboseLevel > 1) {
//      G4cout <<"Momentum direction: (x ,y,z)";
//      G4cout << "(" << momx <<"," << momy << "," << momz << ")" << G4endl;
      }

      // This vector tells the direction of the particle.
      G4ThreeVector momDirection(momx, momy, momz);
      momDirection = momDirection.unit();
        
      // Identify the particle/nucleus:
      G4int particleIdentified = 0;

      // Proton
      if((varntp->avv[particleI] == 1) && (varntp->zvv[particleI] == 1)) {
        cascadeParticle = 
          new G4DynamicParticle(G4Proton::ProtonDefinition(), momDirection, eKin);
        particleIdentified++;
      }

      // Neutron
      if((varntp->avv[particleI] == 1) && (varntp->zvv[particleI] == 0)) {
        cascadeParticle = 
          new G4DynamicParticle(G4Neutron::NeutronDefinition(), momDirection, eKin);
        particleIdentified++;
      }

      // PionPlus
      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == 1)) {
        cascadeParticle = 
          new G4DynamicParticle(G4PionPlus::PionPlusDefinition(), momDirection, eKin);
        particleIdentified++;
      }

      // PionZero
      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == 0)) {
        cascadeParticle = 
          new G4DynamicParticle(G4PionZero::PionZeroDefinition(), momDirection, eKin);
        particleIdentified++;
      }

      // PionMinus
      if((varntp->avv[particleI] == -1) && (varntp->zvv[particleI] == -1)) {
        cascadeParticle = 
          new G4DynamicParticle(G4PionMinus::PionMinusDefinition(), momDirection, eKin);
        particleIdentified++;
      }

      // Nuclei fragment
      if((varntp->avv[particleI] > 1) && (varntp->zvv[particleI] >= 1)) {
        G4ParticleDefinition * aIonDef = 0;
        G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();

        G4int A = G4int(varntp->avv[particleI]);
        G4int Z = G4int(varntp->zvv[particleI]);
        // aIonDef = theTableOfParticles->FindIon(Z, A, 0, Z);
	aIonDef = theTableOfParticles->GetIon(Z, A, varntp->exini*MeV);

	if(aIonDef != 0) {
	  cascadeParticle = 
	    new G4DynamicParticle(aIonDef, momDirection, eKin);
	  particleIdentified++;
	}
      }

      // Check that the particle was identified properly.
      if(particleIdentified == 1) {
        // Put data into G4HadFinalState:
        cascadeParticle->Set4Momentum(cascadeParticle->Get4Momentum()*=toLabFrame);
        theResult.AddSecondary(cascadeParticle); 
      }
      // Particle identification failed. Checking why...
      else {
        // Particle was identified as more than one particle type. 
        if(particleIdentified > 1) {
          G4cout <<"G4InclLightIonInterface: One outcoming particle was identified as";
          G4cout <<"more than one particle type. This is probably due to a bug in the interface." << G4endl;
          G4cout <<"Particle A:" << varntp->avv[particleI] << "Z: " << varntp->zvv[particleI] << G4endl;
          G4cout << "(particleIdentified =" << particleIdentified << ")"  << G4endl;
        }
      }
    }

    // End of conversion

    // Clean up: Clean up the number of generated particles in the
    // common block VARNTP_ for the processing of the next event.
    varntp->ntrack = 0;
    // End of cleanup.
  }
  else {
    // If the bullet type was not recognized by the interface it will
    // be returned back without any interaction.

    theResult.SetStatusChange(stopAndKill);

    G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();
    cascadeParticle = new G4DynamicParticle(theTableOfParticles->FindParticle(aTrack.GetDefinition()), aTrack.Get4Momentum());

    theResult.AddSecondary(cascadeParticle);

    // Unsupported bullets
    G4cout <<"G4InclLightIonInterface: Error processing event number (internal) " << eventNumber << G4endl;

    if(verboseLevel > 3) {
      diagdata <<"G4InclLightIonInterface: Error processing event number (internal) " << eventNumber << G4endl;
    }

    if(bulletType == 0) {
      G4cout <<"G4Incl4AblaCascadeInterface: Unknown bullet type" << G4endl;      
      G4cout <<"Bullet particle name: " << cascadeParticle->GetDefinition()->GetParticleName() << G4endl;
      
      if(verboseLevel > 3) {
        diagdata <<"G4Incl4AblaCascadeInterface: Unknown bullet type" << G4endl;      
        diagdata <<"Bullet particle name: " << cascadeParticle->GetDefinition()->GetParticleName() << G4endl;
      }
    }

    // Unsupported targets
    if((targetA == 1) && (targetZ == 1)) {
      G4cout <<"Unsupported target: " << G4endl;
      G4cout <<"Target A: " << targetA << G4endl;
      G4cout <<"TargetZ: " << targetZ << G4endl;

      if(verboseLevel > 3) {
        diagdata <<"Unsupported target: " << G4endl;
        diagdata <<"Target A: " << targetA << G4endl;
       diagdata <<"TargetZ: " << targetZ << G4endl;
      }
    }

    // Unsupported bullet energies
    if(bulletE < 100) {
      G4cout <<"Unsupported bullet energy: " << bulletE << " MeV. (Lower limit: 100 MeV)." << G4endl;

      if(verboseLevel > 3) {
        diagdata <<"Unsupported bullet energy: " << bulletE << " MeV. (Lower limit: 100 MeV)." << G4endl;
      }
    }
  
    G4cout <<"Failover: returning the original bullet with original energy back to Geant4." << G4endl;
    G4cout <<"G4InclLightIonInterface: End of Error message." << G4endl;

    if(verboseLevel > 3) {
      diagdata <<"Failover: returning the original bullet with original energy back to Geant4." << G4endl;
      diagdata <<"G4InclLightIonInterface: End of Error message." << G4endl;
    }
  }

  // Free allocated memory
  delete varntp;
  delete calincl;
  delete ws;
  delete incl;
  
  return &theResult;
} 

G4ReactionProductVector* G4InclLightIonInterface::Propagate(G4KineticTrackVector* , 
                                                       G4V3DNucleus* ) {
  return 0;
}

