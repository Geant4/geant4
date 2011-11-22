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
#include "G4BinaryCascade.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"
#include "G4StateManager.hh"

#include "G4LeptonConstructor.hh"
#include "G4Gamma.hh"

#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"

#define APPLYYOURSELF 1

/*
#include "ResourceMgr.hh"
// in G4BinaryCascade.cc
ResourceMgr absorbRes;
ResourceMgr buildTargetListRes;
ResourceMgr findCollisionRes;
ResourceMgr findFirstCollisionRes;
ResourceMgr totalRes;
ResourceMgr thePropagatorRes;
ResourceMgr updateRes;
ResourceMgr scattererRes;
void PrintResourceUsage(ResourceMgr & resMgr);
*/
G4V3DNucleus * global3DNucleus;


void Init(G4int code, std::vector<G4ParticleDefinition *> & particles);
void GetParam(G4int & A, G4int & Z, G4double & pMin, G4double & pMax,
	      G4int &nEvents, G4int & particleCode);
G4ThreeVector GetSpherePoint(G4double r);
G4ThreeVector GetSurfacePoint(G4double r);

int main(int argc, char * argv[])
{
    G4StateManager* g4State=G4StateManager::GetStateManager();
    if (! g4State->SetNewState(G4State_Init)) 
      G4cout << "error changing G4state"<< G4endl;;   

  G4String a=G4Proton::ProtonDefinition()->GetParticleName();
  G4cout << " At Start " << a << G4endl;
  G4cout.setf( std::ios::scientific, std::ios::floatfield );
//  totalRes.Start();
  G4LeptonConstructor leptons;
  leptons.ConstructParticle();
  G4IonConstructor ions;
  ions.ConstructParticle();
  const G4ParticleDefinition* proton = G4Proton::Proton();
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  const G4ParticleDefinition* pin = G4PionMinus::PionMinus();
  const G4ParticleDefinition* pip = G4PionPlus::PionPlus();
  const G4ParticleDefinition* pi0 = G4PionZero::PionZero();
  const G4ParticleDefinition* deu = G4Deuteron::DeuteronDefinition();
  const G4ParticleDefinition* tri = G4Triton::TritonDefinition();
  const G4ParticleDefinition* alp = G4Alpha::AlphaDefinition();
  const G4ParticleDefinition* ion = G4GenericIon::GenericIon();
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

  G4int A, Z, nEvents, particleCode;
  G4double pMin, pMax;
  std::vector<G4ParticleDefinition *> particles;
  if(argc == 7)
  {
    A = atoi(argv[1]);
    Z = atoi(argv[2]);
    pMin = G4double(atoi(argv[3]))*MeV;
    pMax = G4double(atoi(argv[4]))*MeV;
    nEvents = atoi(argv[5]);
    particleCode = atoi(argv[6]);
  }
  else
    GetParam(A, Z, pMin, pMax, nEvents, particleCode);

  if(pMin > pMax)
  {
    G4double pTmp = pMin;
    pMin = pMax;
    pMax = pTmp;
  }

  Init(particleCode, particles);

  G4BinaryCascade bic;
  G4ExcitationHandler * excite = new G4ExcitationHandler();
  G4VPreCompoundModel * precompound = new G4PreCompoundModel(excite);

//  bic.SetDeExcitation(0);
  bic.SetDeExcitation(precompound);


  for(G4int event = 0; event < nEvents; ++event)
  {
// Get particle definition
    G4cout << " ---- event " << event << "----------------------"<< G4endl;
    G4double type = G4UniformRand();
    unsigned int index = G4int(type*particles.size());
    if(index == particles.size())  // if random gave 1!
      --index;
    G4ParticleDefinition * particle = particles[index];
// position
    G4ThreeVector pos(0, 0, -5.*fermi);
// momentum
    G4double p = G4UniformRand()*(pMax-pMin)+pMin;
    G4double mass = particle->GetPDGMass();
    G4double e = std::sqrt(mass*mass+p*p);
    G4LorentzVector mom(0, p, 0, e);
    
    G4cout << "initial p /E . Ekin, nucleus " << p << " / " << e 
           << " " << e-mass << " "
	   <<G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A)<< G4endl;

// build the nucleus
#ifdef APPLYYOURSELF
  //  G4Nucleus nucleus(G4double(A), G4double(Z));  // for ApplyYourself()
  G4Nucleus nucleus(A, Z);  // for ApplyYourself()
//  global3DNucleus = new G4Fancy3DNucleus;
//  global3DNucleus->Init(A, Z);
#else
  G4V3DNucleus * fancyNucleus = new G4Fancy3DNucleus;  // for PropagateSTL()
  fancyNucleus->Init(A, Z);
#endif

// build the projectile
#ifdef APPLYYOURSELF
    G4DynamicParticle * dynamic = new G4DynamicParticle(particle, mom);
    G4double time = 0.;
    G4Track projectile(dynamic, time, pos);
    // needed by G4ParticleChange::Initialize(G4Track &)
   // G4Step * step = new G4Step();
    //step->SetStepLength(1.0*mm);
   // projectile.SetStep(step);
    G4HadFinalState * aFinalState = bic.ApplyYourself(projectile, nucleus);

    G4double QValue = 0, Etot=0;
    G4LorentzVector pTot = 0;
//    G4Track * second;
    G4HadSecondary * second;
    const G4DynamicParticle * aSec;
    
    for(G4int isec=0;isec<aFinalState->GetNumberOfSecondaries();isec++)
    {
      second = aFinalState->GetSecondary(isec);
      aSec = second->GetParticle();
      G4cout << "SECONDARIES info ";
      G4cout << aSec->GetDefinition()->GetBaryonNumber() << " ";
      G4cout << aSec->GetTotalEnergy();
      G4cout << aSec->GetMomentum()<< " ";
      G4cout << (1-isec)*aFinalState->GetNumberOfSecondaries();
      G4cout << G4endl;
      QValue += aSec->GetKineticEnergy();
      pTot += aSec->Get4Momentum();
      if ( aSec->GetDefinition()->GetBaryonNumber() == 0 ) 
      {
      	 Etot += aSec->GetTotalEnergy();
      } else {
	 Etot += aSec->GetKineticEnergy();
	 if (aSec->Get4Momentum().mag() > 0.95*GeV )
	 {
	      Etot += aSec->Get4Momentum().mag() - 939*MeV;
	 }
	    
      } 
      

//      delete second;
    }
    G4cout << "QVALUE = "<<QValue<<G4endl;
    G4cout << "Ptotal, pseudoE, exciteE = " 
             << pTot << " " << Etot << " " << e-mass-Etot <<G4endl;

    aFinalState->Clear();
//    delete dynamic;
//    delete step;

//    change->DumpInfo();

#else
    G4double radius = fancyNucleus->GetOuterRadius();
    G4KineticTrackVectorSTL * ktv = new G4KineticTrackVectorSTL;
//    bic.GetSpherePoint(radius, pos); // start from the surface
    pos = GetSpherePoint(radius); // start already inside, pos.mag()< r
    pos = GetSurfacePoint(radius*1.01); // start outside,  pos.mag()= r*1.01
//    pos = G4ThreeVector(1.*fermi,0.*fermi,-15.*fermi);
    G4KineticTrack * kt = new G4KineticTrack(particle, 0., pos, mom);
    ktv->push_back(kt);
    G4ReactionProductVector *result=bic.PropagateSTL(ktv, fancyNucleus);

//    bic.GetExcitationEnergy();
//    std::cout << bic.GetExcitationEnergy()/MeV << endl;
    
    G4ReactionProductVector::iterator result_iter;
    G4int isec=0;
    for(result_iter = result->begin();result_iter != result->end();
        ++result_iter)
    {    
       G4ThreeVector pFinal= (*result_iter)->GetMomentum();

       G4cout << " Final particle # "<< ++isec << ": "
    	     << pFinal.x() << "  "
	     << pFinal.y() << "  "
	     << pFinal.z() << "  "
	     << (*result_iter)->GetTotalEnergy() << "  "
	     << G4endl;
     }	   
#endif
  }
/*
  totalRes.Stop();
  std::cout << "Absorb() resources:" << endl;
  PrintResourceUsage(absorbRes);
  std::cout << "FindFirstCollision() resources:" << endl;
  PrintResourceUsage(findFirstCollisionRes);
  std::cout << "FindCollision() resources:" << endl;
  PrintResourceUsage(findCollisionRes);
  std::cout << "BuildTargetList() resources:" << endl;
  PrintResourceUsage(buildTargetListRes);
  std::cout << "thePropagator() resources:" << endl;
  PrintResourceUsage(thePropagatorRes);
  std::cout << "UpdateTrackAndCollisions() resources:" << endl;
  PrintResourceUsage(updateRes);
  std::cout << "scatterer resources:" << endl;
  PrintResourceUsage(scattererRes);
  std::cout << "main() resources:" << endl;
  PrintResourceUsage(totalRes);
*/
//  std::cout << "Exiting" << endl;
  std::cout.flush();

  return 0;
}

/*

void PrintResourceUsage(ResourceMgr & resMgr)
{
  double userTime = resMgr.GetUserTime();
  double sysTime = resMgr.GetSysTime();
  double cpuTime = userTime+sysTime;
  std::cout << "cpu time: " << cpuTime << " (" << userTime << "+"
       << sysTime << ")" << endl;
}


*/

void Init(G4int code, std::vector<G4ParticleDefinition *> & particles)
{
  if(code == 0)
  {
    particles.push_back(G4Gamma::Gamma());
    return;
  }

  G4int count = 0;
  while(code != 0)
  {
    if(code & 1)
    {
      switch (count)
      {
      case 0:
	particles.push_back(G4Proton::Proton());
	break;
      case 1:
	particles.push_back(G4Neutron::Neutron());
	break;
      case 2:
	particles.push_back(G4AntiProton::AntiProton());
	break;
      case 3:
	particles.push_back(G4PionPlus::PionPlus());
	break;
      case 4:
	particles.push_back(G4PionZero::PionZero());
	break;
      case 5:
	particles.push_back(G4PionMinus::PionMinus());
	break;
      case 6:
	particles.push_back(G4KaonPlus::KaonPlus());
	break;
      case 7:
	particles.push_back(G4KaonZero::KaonZero());
	break;
      case 8:
	particles.push_back(G4KaonMinus::KaonMinus());
	break;
      case 9:
	particles.push_back(G4SigmaPlus::SigmaPlus());
	break;
      case 10:
	particles.push_back(G4SigmaZero::SigmaZero());
	break;
      case 11:
	particles.push_back(G4SigmaMinus::SigmaMinus());
	break;
      }
    }
    ++count;
    code = code>>1;
  }
}


void GetParam(G4int & A, G4int & Z, G4double & pMin, G4double & pMax,
	      G4int &nEvents, G4int & particleCode)
{
// A and Z
  std::cout << "Input Nucleus A and Z: ";
  std::cin >> A >> Z;
// initial impulse
  std::cout << "Input min and max of |momentum|: ";
  std::cin >> pMin >> pMax;
// number of events
  std::cout << "Input number of events: ";
  std::cin >> nEvents;
// particle code
  std::cout << "Binary code for projectile type:" << std::endl;
  std::cout << "1 = proton, 2 = neutron, 4 = antiproton" << std::endl;
  std::cout << "8 = pi+, 16 = pi0, 32 = pi-" << std::endl;
  std::cout << "64 = k+, 128 = k0, 256 = k-" << std::endl;
  std::cout << "512 = s+, 1024 = s0, 2048 = s-" << std::endl;
  std::cout << "       0 for photon" << std::endl;
  std::cout << "Input projectile type (binary code): ";
  std::cin >> particleCode;
}




G4ThreeVector GetSpherePoint(G4double r)
{
// random point INSIDE a sphere

  G4double x = r*(2*G4UniformRand()-1);
  G4double tmp = std::sqrt(r*r-x*x);
  G4double y = tmp*(2*G4UniformRand()-1);
  tmp = std::sqrt(r*r-x*x-y*y);
  G4double z = tmp*(2*G4UniformRand()-1);
  G4ThreeVector point;
  point.setX(x);
  point.setY(y);
  point.setZ(z);
  return point;
}
G4ThreeVector GetSurfacePoint(G4double r)
{
// random point oin surface of a sphere

  G4double x = r*(2*G4UniformRand()-1);
  G4double tmp = std::sqrt(r*r-x*x);
  G4double y = tmp*(2*G4UniformRand()-1);
  G4double z = -std::sqrt(r*r-x*x-y*y);
  G4ThreeVector point;
  point.setX(x);
  point.setY(y);
  point.setZ(z);
  return point;
}
