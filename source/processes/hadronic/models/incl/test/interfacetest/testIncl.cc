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
// $Id: testIncl.cc,v 1.4 2007-10-31 14:57:24 gcosmo Exp $

#include <vector>

#include <iomanip>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

#include "globals.hh"
#include "Randomize.hh"

#include "G4StateManager.hh"
#include "G4ParticleTable.hh"

#include "G4Collider.hh"
#include "G4InuclCollider.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4CollisionOutput.hh"
#include "G4Analyser.hh"
#include "G4WatcherGun.hh"
#include "G4ios.hh"
#include "G4BertiniData.hh"
#include "G4CascadSpecialFunctions.hh"
#include "G4IonTable.hh"
#include "G4Nucleus.hh"
#include "G4NucleiModel.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4Proton.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4HadFinalState.hh"
#include "G4DynamicParticle.hh"

#include "G4CascadeInterface.hh"
#include "G4ElasticCascadeInterface.hh"
//#include "G4InuclEvaporation.hh"

#include "G4InclCascadeInterface.hh"
#include "G4InclAblaCascadeInterface.hh"

#include "G4RunManager.hh"

void test(std::string, int);

G4int tCascadeInterface();

int main(int argc, char **argv ) {
  G4int verboseLevel = 2;
  G4cout << "Geant4 cascade region benchmarks" << G4endl;

  if (argc < 2)
    {
      printf("usage: benchmarks <test ID> <parameters>\n");
      return(1);
    }
  enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

  // defaults
  G4int testId      = 1;
  G4int nCollisions = 10;      // collisions to be generated
  G4int  bulletType = proton;    // bullet particle
  G4double     momZ = 160;      // momentum in z-direction
  G4double        A = 27.0;      // target atomic weight Al
  G4double        Z = 13.0;      // target atomic number
  testId      =           (argc > 2) ? atoi(argv[1]) : testId;
  nCollisions =           (argc > 2) ? atoi(argv[1]) : nCollisions;
  bulletType  =           (argc > 3) ? atoi(argv[2]) : proton;
  momZ        = G4double(((argc > 4) ? atoi(argv[3]) : momZ)) / GeV;
  A           = G4double(((argc > 5) ? atoi(argv[4]) : A));
  Z           = G4double(((argc > 7) ? atoi(argv[5]) : Z));

  if (verboseLevel > 1) {
    G4cout << " nCollisions " << nCollisions << G4endl;
    G4cout << "  bulletType " << bulletType  << G4endl;
    G4cout << "        momZ " << momZ        << G4endl;
    G4cout << "           A " << A           << G4endl;
    G4cout << "           Z " << Z           << G4endl;
  }

  test("Cascade interface", tCascadeInterface());

  return 0;       
}

void test(std::string txt, int testStatus) {
  G4cout << txt << ": ";
  if (testStatus){ 
    G4cout << "OK";
  } else {
    G4cout << "Fail" << G4endl;
  }; 

  G4cout << G4endl;  // test timing 
}


int tCascadeInterface() {
  G4int verboseLevel                = 2;                          
  if (verboseLevel > 1) {
    G4cout << ">>> tCascadeInterface start" << G4endl;
  }
  G4int numberOfCascades            = 1; 
  G4double projectileMomentum       = 1.0 * GeV;

  //G4ParticleDefinition *particle = G4PionMinus::PionMinus();  
  //G4ParticleDefinition *particle = G4Neutron::Neutron();  
  G4ParticleDefinition *particle = G4Proton::Proton();  

  G4DynamicParticle *projectile = new G4DynamicParticle(); 
  projectile->SetDefinition(particle);
  //projectile->SetKineticEnergy( 1.0 * GeV);
  projectile->SetMomentum(projectileMomentum);
  projectile->SetMomentumDirection(1.0, 0.0, 0.0);  

  if (verboseLevel > 1) {
    G4cout << "projectile" << G4endl;
    G4cout << " type           : " << projectile->GetDefinition()        << G4endl;
    G4cout << " kinetic energy : " << projectile->GetKineticEnergy()     << G4endl;
    G4cout << " momentum       : " << projectile->GetMomentum()          << G4endl;
    G4cout << " p direction    : " << projectile->GetMomentumDirection() << G4endl;
  }

  // Set projectile particle track
  G4ThreeVector v;                                            
  v.setX(0.0 * fermi); 
  v.setY(0.0 * fermi); 
  v.setZ(0.0 * fermi);
  G4Track aTrack(projectile, 0, v);

  // Set target nucleus
  G4Nucleus targetNucleus;                                        
//   G4double a(10);
//   G4double z(10);
   G4double a(207);
   G4double z(82);
   //   G4double a(27);
   //   G4double z(13);

  targetNucleus.SetParameters(a, z);

  if (verboseLevel > 1) {
    G4cout << "target" << G4endl;
    G4cout << " a              : " << a                              << G4endl;
    G4cout << " z              : " << z                              << G4endl;
    G4cout << " atomic mass    : " << targetNucleus.AtomicMass(a, z) << G4endl;
  }

  //  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& theNucleus); 
  G4double            px = 0;
  G4double            py = 0;
  //  G4double            pz = 2000;
  G4double            pz = 2000;
  G4ThreeVector       inVector(px, py, pz);
  G4ThreeVector       outVector, aPosition(0., 0., 0.);

  G4Proton          * aProton   = G4Proton::Proton();

  G4DynamicParticle   aParticle;
  aParticle.SetDefinition(aProton);
  //  G4double Momentum = 2000;
  G4double Momentum = 2000;
  G4double Tkin = std::sqrt(Momentum*Momentum+938.27*938.27)-938.27;

  inVector.setZ(Momentum);

  G4cout <<  "  inVector " <<inVector.x()<<" " << 
    inVector.y()<<" "<<  inVector.z()<<" " << " Tkin " << Tkin <<G4endl;

  aParticle.SetMomentum(inVector);

  aParticle.SetKineticEnergy(Tkin);

  const  G4HadProjectile hadProj = G4HadProjectile(aParticle);

  G4HadFinalState   * hadSta    = new G4HadFinalState();   

  // Set the particle table singleton ready:
  //  G4ParticleTable::GetParticleTable()->SetReadiness();

  //G4CascadeInterface *theCascade  = new G4CascadeInterface();
  G4InclCascadeInterface *theCascade = new G4InclCascadeInterface();
 
  for (G4int cascadeID =1 ; cascadeID <= numberOfCascades; cascadeID++) { 
    if (verboseLevel > 1) G4cout << "inc " << cascadeID << G4endl;
    hadSta = theCascade->ApplyYourself(hadProj, targetNucleus);

    G4int nPart = hadSta->GetNumberOfSecondaries();
    G4cout << "  # secondaries " << nPart << G4endl;
    outVector  =  hadSta->GetMomentumChange();
    G4cout << "  momentum change " << outVector << G4endl;
    G4double outE = hadSta->GetEnergyChange();
    G4cout << "  energy change " << outE << G4endl;
    G4double outP = std::sqrt(outE*outE-938.27*938.27);

    for (G4int iSecondary =1 ; iSecondary < nPart; iSecondary++) { 

      G4HadSecondary    * NuclSecond = hadSta->GetSecondary(iSecondary);
      G4cout << "    secondary         " << iSecondary << G4endl;  
      G4DynamicParticle * secPart = NuclSecond->GetParticle();

      G4cout<<"      nucleus name      " << secPart->GetDefinition()->GetParticleName() << G4endl;

      G4ThreeVector outVectorN =  secPart->GetMomentum();
      G4cout << "      out vector      " << outVectorN << G4endl;
      G4double outEtot =  secPart->GetTotalEnergy();
      G4cout << "      particle  tot E " << outEtot << G4endl;  
      G4double outEkin =  secPart->GetKineticEnergy();
      G4cout << "      particle  kin E " << outEkin << G4endl;  
    }
  }

  delete projectile;
  delete theCascade;
  return 1;   
}
