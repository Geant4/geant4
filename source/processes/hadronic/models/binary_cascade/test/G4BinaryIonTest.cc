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
#include "G4BinaryLightIonReaction.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"

#include "G4LeptonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Gamma.hh"

G4V3DNucleus * global3DNucleus;


void Init(G4int code, vector<G4ParticleDefinition *> & particles);
void GetParam(G4int & A, G4int & Z, G4double & pMin, G4double & pMax,
	      G4int &nEvents, G4int & particleCode);
G4ThreeVector GetSpherePoint(G4double r);
G4ThreeVector GetSurfacePoint(G4double r);

int main(int argc, char * argv[])
{
  G4String a=G4Proton::ProtonDefinition()->GetParticleName();
  G4cout << " At Start " << a << G4endl;
  G4cout.setf( std::ios::scientific, std::ios::floatfield );
//  totalRes.Start();
  G4LeptonConstructor leptons;
  leptons.ConstructParticle();
  G4IonConstructor ions;
  ions.ConstructParticle();

  G4int A, Z, nEvents, particleCode;
  G4double pMin, pMax;
  vector<G4ParticleDefinition *> particles;
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

  G4BinaryCascade hkm;
  hkm.SetDeExcitation(0);

  for(G4int event = 0; event < nEvents; ++event)
  {
    G4cerr << " ---- event " << event << "----------------------"<< G4endl;
    G4V3DNucleus * fancyNucleus = new G4Fancy3DNucleus;  // for Propagate()
    fancyNucleus->Init(16, 8);
    G4V3DNucleus * projectile = new G4Fancy3DNucleus;
    projectile->Init(12, 6);
    G4double radius = fancyNucleus->GetOuterRadius();
    G4double impactMax = fancyNucleus->GetOuterRadius()+projectile->GetOuterRadius();
// position
    G4double aX=(2.*G4UniformRand()-1.)*impactMax;
    G4double aY=(2.*G4UniformRand()-1.)*impactMax;
    
    G4ThreeVector pos(aX, aY, -20.*fermi);
// momentum
    G4double mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(6, 12);
    G4double p = 280*MeV;
    G4double e = std::sqrt(mass*mass+p*p);
    G4LorentzVector mom(0, p, 0, e);
    
// build the projectile
    G4KineticTrackVector * ktv = new G4KineticTrackVector;
    G4ThreeVector beta = mom.beta();
    projectile->DoLorentzBoost(beta);
    projectile->StartLoop();
    G4Nucleon * aNuc;
    while( (aNuc=projectile->GetNextNucleon()) )
    {
      G4KineticTrack * it = 
          new G4KineticTrack(aNuc->GetDefinition(), 0, 
	                     aNuc->GetPosition()+pos, const_cast<G4LorentzVector &>(aNuc->GetMomentum()) );
      ktv->push_back(it);
    }

// do the reactions
    G4ReactionProductVector *result=hkm.Propagate(ktv, fancyNucleus);
    
    delete fancyNucleus;
    delete projectile;
    
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
  }
  cout.flush();
  return 0;
}

/*

void PrintResourceUsage(ResourceMgr & resMgr)
{
  double userTime = resMgr.GetUserTime();
  double sysTime = resMgr.GetSysTime();
  double cpuTime = userTime+sysTime;
  cout << "cpu time: " << cpuTime << " (" << userTime << "+"
       << sysTime << ")" << endl;
}


*/

void Init(G4int code, vector<G4ParticleDefinition *> & particles)
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
  cout << "Input Nucleus A and Z: ";
  cin >> A >> Z;
// number of events
  cout << "Input number of events: ";
  cin >> nEvents;
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
