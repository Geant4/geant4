// J.P. Wellisch, X-mas 2002.

#include "G4Rutherford.hh"
#include "G4HadPointVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4LorentzRotation.hh"

G4VParticleChange * G4Rutherford::
ApplyYourself(const G4Track& aTrack, G4Nucleus& targetNucleus)
{
  theResult.Initialize(aTrack);
  theResult.SetStatusChange(fStopAndKill);
  theResult.SetNumberOfSecondaries(2);
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4double m1=aParticle->GetDefinition()->GetPDGMass();
  G4double z2=targetNucleus.GetZ();
  G4double a2=targetNucleus.GetN();
  G4double m2=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z2,a2);

  G4LorentzVector projectileMomentum = aTrack.GetDynamicParticle()->Get4Momentum();
  G4LorentzRotation toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4LorentzRotation toLabFrame = toZ.inverse();

  if(game && m1)
  {
    // all integrals prepared.
  }
  else
  {
    G4double cth = Apply(aParticle->GetDefinition(), targetNucleus); // cos(th) of the lighter
    G4double th = acos(cth);
    G4double phi = 2.*pi*G4UniformRand();
    if(m2>m1)
    {
      // get momentum
      G4double p = projectileMomentum.vect().mag();

      G4double a = m1/m2+2.*pow(sin(th), 2.);
      G4double b = 2.*p*cth;
      G4double c =p*p*(1-m1/m2);
      G4double p1=(-b+sqrt(b*b-4.*a*c))/(2.*a);
      // transform to lab.
      G4ThreeVector it(p1*sin(th)*sin(phi), p1*sin(th)*cos(phi), p1*cos(th));
      G4LorentzVector theMom(it, sqrt(m1*m1+it.mag2()) );
      theMom*=toLabFrame;
      it=theMom.vect();
      G4DynamicParticle * aSec = 
	  new G4DynamicParticle(aParticle->GetDefinition(), it.unit(), it.mag2()/(2.*m1));
      theResult.AddSecondary(aSec);
      G4ThreeVector it1(-it.x(), -it.y(), p-it.z());
      G4LorentzVector theMom1(it1, sqrt(m2*m2+it1.mag2()) );
      theMom1*=toLabFrame;
      it1 = theMom1.vect();
      G4ParticleDefinition * theSec = 
          G4ParticleTable::GetParticleTable()->FindIon(z2, a2, 0, z2);
      G4DynamicParticle * bSec = 
	  new G4DynamicParticle(theSec, it1.unit(), it1.mag2()/(2.*m2));
      theResult.AddSecondary(bSec);
    }
    else
    {
      // get momentum
      G4double p = projectileMomentum.vect().mag();
      G4double sth=sin(th);

      G4double p2=-2.*p*cth;
      p2/=(cth*cth-sth*sth+m1/m2);
      // transform to lab.
      G4double vPrim = projectileMomentum.vect().mag()/m1;
      G4ThreeVector it(p2*sin(th)*sin(phi), p2*sin(th)*cos(phi), p2*cos(pi-th));
      it.setZ(it.z()+vPrim*m2);
      G4LorentzVector theMom(it, sqrt(m2*m2+it.mag2()) );
      theMom*=toLabFrame;
      it=theMom.vect();
      G4ParticleDefinition * theSec = 
          G4ParticleTable::GetParticleTable()->FindIon(z2, a2, 0, z2);
      G4DynamicParticle * aSec = 
	  new G4DynamicParticle(theSec, it.unit(), it.mag2()/(2.*m2));
      theResult.AddSecondary(aSec);
      G4ThreeVector it1(-it.x(), -it.y(), p-it.z());
      G4LorentzVector theMom1(it1, sqrt(m1*m1+it1.mag2()) );
      theMom1*=toLabFrame;
      it1 = theMom1.vect();
      G4DynamicParticle * bSec = 
	  new G4DynamicParticle(aParticle->GetDefinition(), it1.unit(), it1.mag2()/(2.*m1));
      theResult.AddSecondary(bSec);
    }
  }
  return &theResult;
}

G4double G4Rutherford::
Apply(const G4ParticleDefinition* aP, G4Nucleus& targetNucleus)
{
  G4double m1=aP->GetPDGMass();
  G4double z2=targetNucleus.GetZ();
  G4double a2=targetNucleus.GetN();
  G4double m2=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z2,a2);
  vector<G4double> theRuther;
  vector<G4double> theValue;
  G4double start = theMaxCosTh;
  G4double eps = 0.0001;
  G4double r=m1/m2;
  if(r>1) r=1./r;
  G4double ruth=0;
  G4int count=0;
  while(start>-1)
  {
    ruth = Ruther(r, start);
    theRuther.push_back(ruth);
    theValue.push_back(start);
    count ++;
    start -=eps;
    eps *= 1.03;
  }
  ruth = Ruther(r, -.99999999);
  theRuther.push_back(ruth);
  theValue.push_back(-1);

  vector<G4double> integral;
  double inte = 0;
  integral.push_back(0);
  size_t i=0;
  for(i=0; i<theRuther.size()-1; i++)
  {
    double y1=theRuther[i];
    double y2=theRuther[i+1];
    double x1=theValue[i];
    double x2=theValue[i+1];
    inte += 0.5*(y1+y2)*(x2-x1);
    integral.push_back(inte);
//    cout << theValue[i+1] <<" "<<integral[i+1]<<endl;
  }
//  cout << integral.size()<<" "<<theRuther.size()<<endl;
  for(i=0; i<integral.size(); i++)
  {
    integral[i]/=inte;
//    cout << theValue[i] <<" "<<integral[i]<<endl;
  }
  // fill particle change
  G4double cth=0;
//  for (i=0; i<1000000; i++)
//  {
    G4double random=G4UniformRand();
    size_t j;
    for(j=0; j<integral.size(); j++)
    {
      if(integral[j]>random) break;
    }
    double x1=theValue[j-1];
    double x2=theValue[j];
    double y1=integral[j-1];
    double y2=integral[j];
    cth = x1+(random-y1)*(x2-x1)/(y2-y1);
    G4cout << x1 << G4endl;
//  }
  return cth;
}
