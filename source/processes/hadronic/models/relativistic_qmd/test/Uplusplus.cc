
#include "G4KineticTrackVector.hh"
#include "G4KineticTrack.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4UppModel.hh"
#include "G4UppSimpleFieldtransport.hh"
#include "G4VUppAnalyzer.hh"
#include "G4Proton.hh"
#include <iostream>


class SimpleAnalyzer : public G4VUppAnalyzer
{
public:
  SimpleAnalyzer(const string& ofile) : outputfile(ofile) { }
  void analyze(const G4UppTrackVector& all) const;
  void analyze(const G4UppTrackVector& all, const G4UppTrackChange& i) const
    {}
  string getName() const { return "SimpleAnalyzer"; }
private:
  const string outputfile;
};



void SimpleAnalyzer::analyze(const G4UppTrackVector& all) const
{
  string fnam;
  cout << "Analyzer running..." << endl;
  cout << "type filename: ";
  cin >> fnam;
  ofstream out(fnam.c_str());
  for (int i=0; i<all.size(); i++) {
    G4ThreeVector aPos = all[i]->GetPosition();
    out << aPos.x()/fermi << " " << aPos.y()/fermi << " " << aPos.z()/fermi << endl;
  }
}



G4int main()
{
  G4UppModel myQMD;
  SimpleAnalyzer myAnalyzer("step0");
  G4KineticTrackVector* final;
  G4UppSimpleFieldtransport aFieldtransport;

  cout << endl;
  cout << "U++ Version 0.0.2" << endl;
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
  // myQMD.addAnalyzer(&ana, -5.0*fermi);
  // myQMD.addAnalyzer(&ana,  0.0*fermi);
  // myQMD.addAnalyzer(&ana,  5.0*fermi);
  // myQMD.addAnalyzer(&ana, 15.0*fermi);

  G4double ebeam = 3*GeV;
  G4double m0 = G4Proton::ProtonDefinition()->GetPDGMass();
  G4double plab = sqrt(ebeam*(ebeam+2*m0));
  G4double betacm = plab/(ebeam+2*m0);
  cout << "(debug) p_lab: " << plab/GeV << endl;
  cout << "(debug) beta_cm: " << betacm << endl;

  cout << "generating projectile ... " << endl;
  G4Fancy3DNucleus Projectile;
  Projectile.Init(4,2);
  Projectile.DoLorentzBoost(G4ThreeVector(0, 0, betacm));
  Projectile.DoLorentzContraction(G4ThreeVector(0, 0, betacm));
  Projectile.DoTranslation(G4ThreeVector(0, 0, -10.0*fermi));

  cout << "generating target ... " << endl;
  G4Fancy3DNucleus Target;
  Target.Init(4,2);
  Target.DoLorentzBoost(G4ThreeVector(0, 0, -betacm));
  Target.DoLorentzContraction(G4ThreeVector(0, 0, -betacm));
  Target.DoTranslation(G4ThreeVector(0, 0, 10*fermi));

  cout << "initialize model ... " << endl;
  myQMD.initialize(Projectile,Target);

  // myQMD.addAnalyzer(&ana,  0.0*fermi);
  // myQMD.addAnalyzer(&ana,  5.0*fermi);
  // myQMD.addAnalyzer(&ana, 15.0*fermi);

  cout << "starting propagation ... " << endl;
  myQMD.propagate(aFieldtransport);

  // final = myQMD.GetFinalState();

  cout << "finished." << endl;

  return 0;
}
