
#include "G4KineticTrackVector.hh"
#include "G4KineticTrack.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4UppModel.hh"
#include "G4UppSimpleFieldtransport.hh"
#include "G4VUppAnalyzer.hh"
#include "G4Proton.hh"
//#include "G4PionZero.hh"
//#include "G4MuonPlus.hh"
#include <iostream>

#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"

#include "SimpleAnalyzer.hh"



G4int main()
{
  G4UppModel myQMD;

  // G4ShortLivedConstructor ShortLived;
  // G4ShortLivedTable ShortLivedTab;
  // ShortLived.ConstructParticle();
  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();
  G4MesonConstructor Mesons;
  Mesons.ConstructParticle();

  SimpleAnalyzer myAnalyzer("step0");
  G4KineticTrackVector* final;
  G4UppSimpleFieldtransport aFieldtransport;
  //G4ParticleDefinition* aPi0 = G4PionZero::PionZero();
  //G4ParticleDefinition* aMuonPlus = G4MuonPlus::MuonPlus();
  cout << endl;
  cout << "U++ Version 0.0.3" << endl;
  cout << "-----------------" << endl;
  cout << endl;

  // G4KineticTrack* t1 = new G4UppTrack;
  // t1->SetPosition(G4ThreeVector(-1.89332675*fermi,-0.838326437*fermi,-1.19810285*fermi));
  // t1->Set4Momentum(G4LorentzVector(0.0,  0.0,  6.84835755*GeV,6.9094888*GeV) );
  // t1->SetDefinition(G4Proton::Proton());
  // G4KineticTrack* t2 = new G4UppTrack;
  // t2->SetPosition(G4ThreeVector(-3.43756102*fermi ,-2.50154287*fermi , 0.519287767*fermi));
  // t2->Set4Momentum(G4LorentzVector(0.,  0., -6.84835755*GeV, 6.9097837*GeV));
  // t2->SetDefinition(G4Proton::Proton());
  // G4KineticTrackVector t;
  // t.append(t1);
  // t.append(t2);
  // myQMD.Initialize(t);
  // myQMD.addAnalyzer(&ana, -5.0*fermi/c_light);
  // myQMD.addAnalyzer(&ana,  0.0*fermi/c_light);
  // myQMD.addAnalyzer(&ana,  5.0*fermi/c_light);
  // myQMD.addAnalyzer(&ana, 15.0*fermi/c_light);

  G4double ebeam = 2.0*GeV;
  G4double m0 = G4Proton::ProtonDefinition()->GetPDGMass();
  G4double plab = sqrt(ebeam*(ebeam+2*m0));
  G4double betacm = plab/(ebeam+2*m0);
  G4double b_imp = 3*fermi;

  cout << "(debug) p_lab: " << plab/GeV << endl;
  cout << "(debug) beta_cm: " << betacm << endl;

  cout << "generating projectile ... " << endl;
  G4Fancy3DNucleus Projectile;
  // Projectile.Init(1,1);
  Projectile.Init(197,98);
  Projectile.DoLorentzBoost(G4ThreeVector(0, 0, betacm));
  Projectile.DoLorentzContraction(G4ThreeVector(0, 0, betacm));
  Projectile.DoTranslation(G4ThreeVector(b_imp, 0, -10.0*fermi));

  cout << "generating target ... " << endl;
  G4Fancy3DNucleus Target;
  // Target.Init(1,1);
  Target.Init(197,98);
  Target.DoLorentzBoost(G4ThreeVector(0, 0, -betacm));
  Target.DoLorentzContraction(G4ThreeVector(0, 0, -betacm));
  Target.DoTranslation(G4ThreeVector(-b_imp, 0, 10.0*fermi));

  cout << "initialize model ... " << endl;
  myQMD.initialize(Projectile,Target);

  myQMD.addAnalyzer(myAnalyzer, 10.0*fermi/c_light);
  // myQMD.addAnalyzer(myAnalyzer,  5.0*fermi/c_light);
  // myQMD.addAnalyzer(myAnalyzer, 10.0*fermi/c_light);
  myQMD.addAnalyzer(myAnalyzer, 15.0*fermi/c_light);
  myQMD.addAnalyzer(myAnalyzer, 20.0*fermi/c_light);

  cout << "starting propagation ... " << endl;
  myQMD.propagate(aFieldtransport);

  // final = myQMD.GetFinalState();

  cout << "finished." << endl;

  return 0;
}
