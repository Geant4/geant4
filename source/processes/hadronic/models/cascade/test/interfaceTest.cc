//#include <stdio.h>

#include "G4ios.hh"
#include "G4CascadeInterface.hh"
#include "G4ParticleTable.hh"
#include "G4BaryonConstructor.hh" 
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4LeptonConstructor.hh"
//#include "G4Fancy3DNucleus.hh"

//#include "myg4templates.hh"

void ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4BosonConstructor Bosons;
  Bosons.ConstructParticle();

  G4MesonConstructor Mesons;
  Mesons.ConstructParticle();

  G4BaryonConstructor Baryons;
  Baryons.ConstructParticle();

  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  
  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();
}

void KineticTrackVectorInfo(G4KineticTrackVector& KTV, G4LorentzVector* Mom, G4double* Q);

struct
{
  double x;
  double y;
  double z;
  double t;
  double Px;
  double Py;
  double Pz;
  double e;
  double mass;
  int    Encoding;
  int    nHadron;
  int    nEvent;
  int    cHadron;
} Buff;

void main() 
{
  ConstructParticle();
  G4CascadeInterface model;
       
  G4double At = 1.0;
  G4double Zt = 1.0;
  G4double Po = 100 * MeV;
    
  //  G4cout << "Usage:" << G4endl << "KinModelTest A_target, Z_target, P_projectile" << G4endl;

  G4cout << "Enter A_target (1-250) :" << G4endl;
  G4cin >> At;

  G4cout << "Enter Z_target (1-Atarget) :";
  G4cin >> Zt;

  G4cout << "Enter kinetic energy of projectile in MeV:";
  G4cin >> Po;

  G4int nevents;
  G4cout << "Enter number of events:";
  G4cin >> nevents;
   
  if(nevents > 1000) nevents = 1000;
  if(nevents <= 0)   nevents = 10;
   
  char * Ptr = "kyky.dat";

  FILE* out = fopen(Ptr, "wb");
  if(!out) 
    {
      cout << "I cannot open the output file" << G4endl;
      return;
    }

  if(At < 1.0) At = 1.0;
  if(Zt > At) Zt = At;
  if(Po <= 0.0) Po = 100 * MeV; 
  else Po *= MeV;
                                                                                  
  // Create the target
  G4Nucleus theTarget(At, Zt);

  // Create the projectile
  G4ThreeVector position;
  G4ParticleDefinition *pd = G4ParticleTable::GetParticleTable()->FindParticle(2212); // 2212 proton
  // 2112 neutron
  G4double m = pd->GetPDGMass();
  G4double e = m + Po;
  Po = sqrt(e * e - m * m);   
  G4ThreeVector projectileMom(0, 0, Po);
  G4Track theProjectile(new G4DynamicParticle(pd, projectileMom), 0.0 , position);
  G4LorentzVector InMom(theProjectile.GetMomentum(), theProjectile.GetTotalEnergy()
			+ theTarget.AtomicMass(At, Zt));
   
  for(G4int cEvent=1; cEvent<=nevents; cEvent++) {
    G4VParticleChange * theParticleChange = model.ApplyYourself(theProjectile, theTarget);
      
    // Convert output into KineticTrackVector
    G4KineticTrackVector ktvOutput;
    G4int nSecondaries = theParticleChange->GetNumberOfSecondaries();
    for(G4int is=0; is<nSecondaries; is++) {
      G4Track *t = theParticleChange->GetSecondary(is);
      G4ThreeVector pos = t->GetPosition();
      G4LorentzVector mom4(t->GetMomentum(), t->GetTotalEnergy());
      ktvOutput.insert(new G4KineticTrack(
					  t->GetDefinition(),                 // particle definition
					  0.0,                                // formation time
					  pos,                                // position 3-vector
					  mom4));                             // momentum 4-vector
            
      delete t;
    }

    theParticleChange->Clear();
    //delete theParticleChange;

    G4LorentzVector OutMom;
    G4double OutCharge = 0;

    KineticTrackVectorInfo(ktvOutput, &OutMom, &OutCharge);

    cout << "Input 4-momentum" <<  InMom << G4endl;
    cout << "Input sum of charge " << 1+Zt << G4endl;

    cout << "Output 4-momentum" << OutMom << G4endl;
    cout << "Output sum of charge " << OutCharge << G4endl << G4endl;

    G4int nHadron = ktvOutput.length();
    for(G4int c1=0; c1<nHadron; c1++) {
      G4KineticTrack*  KT = ktvOutput.at(c1);
      G4ThreeVector    Pos = KT->GetPosition();
      G4LorentzVector  Mom = KT->Get4Momentum();
      Buff.x  = Pos.x();  Buff.y  = Pos.y();  Buff.z  = Pos.z();
      Buff.t  = KT->GetFormationTime();
      Buff.Px = Mom.px(); Buff.Py = Mom.py(); Buff.Pz = Mom.pz();
      Buff.e  = Mom.e ();
      Buff.mass     = KT->GetDefinition()->GetPDGMass();
      Buff.Encoding = KT->GetDefinition()->GetPDGEncoding();
      Buff.nEvent   = cEvent;
      Buff.nHadron  = nHadron;
      Buff.cHadron  = c1;
      fwrite(&Buff, sizeof(Buff), 1, out);
    }       

    ktvOutput.clearAndDestroy();
  }
   
  fclose(out);
}

void KineticTrackVectorInfo(G4KineticTrackVector& KTV, G4LorentzVector* Mom, G4double* Q)
{
  for(G4int c1 = 0; c1 < KTV.length(); c1++)
    {
      G4KineticTrack& KT =*(KTV.at(c1));
      *Mom += KT.Get4Momentum();
      *Q   += KT.GetDefinition()->GetPDGCharge();
    }
}
