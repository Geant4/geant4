// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPElasticFS.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleTable.hh"
#include "G4NeutronHPDataUsed.hh"

  void G4NeutronHPElasticFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
  {
    G4String tString = "/FS/";
    G4bool dbool;
    G4NeutronHPDataUsed aFile = theNames.GetName(A, Z, dirName, tString, dbool);
    G4String filename = aFile.GetName();
    theBaseA = aFile.GetA();
    theBaseZ = aFile.GetZ();
    if(!dbool)
    {
      hasAnyData = false;
      hasFSData = false; 
      hasXsec = false;
      return;
    }
    G4std::ifstream theData(filename, G4std::ios::in);
    theData >> repFlag >> targetMass >> frameFlag;
    if(repFlag==1)
    {
      G4int nEnergy;
      theData >> nEnergy; 
      theCoefficients = new G4NeutronHPLegendreStore(nEnergy);
      theCoefficients->InitInterpolation(theData);
      G4double temp, energy;
      G4int tempdep, nLegendre;
      G4int i, ii;
      for (i=0; i<nEnergy; i++)
      {
        theData >> temp >> energy >> tempdep >> nLegendre;
        energy *=eV;
        theCoefficients->Init(i, energy, nLegendre);
        theCoefficients->SetTemperature(i, temp);
        G4double coeff=0;
        for(ii=0; ii<nLegendre; ii++)
        {
          // load legendre coefficients.
          theData >> coeff;
          theCoefficients->SetCoeff(i, ii+1, coeff); // @@@HPW@@@
        }
      }
    }
    else if (repFlag==2)
    {
      G4int nEnergy;
      theData >> nEnergy;
      theProbArray = new G4NeutronHPPartial(nEnergy, nEnergy);
      theProbArray->InitInterpolation(theData);
      G4double temp, energy;
      G4int tempdep, nPoints;
      for(G4int i=0; i<nEnergy; i++)
      {
        theData >> temp >> energy >> tempdep >> nPoints;
        energy *= eV;
        theProbArray->InitInterpolation(i, theData);
        theProbArray->SetT(i, temp);
        theProbArray->SetX(i, energy);
        G4double prob, costh;
        for(G4int ii=0; ii<nPoints; ii++)
        {
          // fill probability arrays.
          theData >> costh >> prob;
          theProbArray->SetX(i, ii, costh);
          theProbArray->SetY(i, ii, prob);
        }
      }
    }
    else if (repFlag==0)
    {
      theData >> frameFlag;
    }
    else
    {
      G4cout << "unusable number for repFlag: repFlag="<<repFlag<<G4endl;
      G4Exception("G4NeutronHPElasticFS::Init -- unusable number for repFlag");
    }
  }
  G4ParticleChange * G4NeutronHPElasticFS::ApplyYourself(const G4Track & theTrack)
  {  
    G4int i, ii, iii;
//    G4cout << "G4NeutronHPElasticFS::ApplyYourself+"<<G4endl;
    theResult.Initialize(theTrack);   
    G4double eKinetic = theTrack.GetKineticEnergy();
    const G4DynamicParticle *incidentParticle = theTrack.GetDynamicParticle();
    G4ReactionProduct theNeutron( incidentParticle->GetDefinition() );
    theNeutron.SetMomentum( incidentParticle->GetMomentum() );
    theNeutron.SetKineticEnergy( eKinetic );
//    G4cout << "G4NeutronHPElasticFS::ApplyYourself++"<<eKinetic<<" "<<G4endl;
//    G4cout << "CMSVALUES 0 "<<theNeutron.GetTotalMomentum()<<G4endl;
    G4double pold = theNeutron.GetTotalMomentum();
    
    G4ReactionProduct theTarget; 
    G4Nucleus aNucleus;
    theTarget = aNucleus.GetThermalNucleus( targetMass );
//     G4cout << "Nucleus-test"<<" "<<targetMass<<" ";
//     G4cout << theTarget.GetMomentum().x()<<" ";
//     G4cout << theTarget.GetMomentum().y()<<" ";
//     G4cout << theTarget.GetMomentum().z()<<G4endl;
    
    // neutron and target defined as reaction products.

// prepare lorentz-transformation to Lab.

    G4ThreeVector the3Neutron = theNeutron.GetMomentum();
    G4double nEnergy = theNeutron.GetTotalEnergy();
    G4ThreeVector the3Target = theTarget.GetMomentum();
//    cout << "@@@" << the3Target<<G4endl;
    G4double tEnergy = theTarget.GetTotalEnergy();
    G4ReactionProduct theCMS;
    G4double totE = nEnergy+tEnergy;
    G4ThreeVector the3CMS = the3Target+the3Neutron;
    theCMS.SetMomentum(the3CMS);
    G4double cmsMom = sqrt(the3CMS*the3CMS);
    G4double sqrts = sqrt((totE-cmsMom)*(totE+cmsMom));
    theCMS.SetMass(sqrts);
    theCMS.SetTotalEnergy(totE);
    
    // data come as fcn of n-energy in nuclear rest frame
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutron, theTarget);
    eKinetic = boosted.GetKineticEnergy(); // get kinetic energy for scattering
    G4double cosTh = -2;
    if(repFlag == 1)
    {
      cosTh = theCoefficients->SampleElastic(eKinetic);
    }
    
    else if (repFlag==2)
    {
      cosTh = theProbArray->Sample(eKinetic);
    }
    else if (repFlag==0)
    {
      cosTh = 2.*G4UniformRand()-1.;
    }
    else
    {
      G4cout << "unusable number for repFlag: repFlag="<<repFlag<<G4endl;
      G4Exception("G4NeutronHPElasticFS::Init -- unusable number for repFlag");
    }
    if(cosTh<-1.1) return NULL;
    G4double phi = twopi*G4UniformRand();
    G4double theta = acos(cosTh);
    G4double sinth = sin(theta);
    if (frameFlag == 1) // final state data given in target rest frame.
    {
      // we have the scattering angle, now we need the energy, then do the
      // boosting.
      // relativistic elastic scattering energy angular correlation:
      theNeutron.Lorentz(theNeutron, theTarget);
      G4double e0 = theNeutron.GetTotalEnergy();
      G4double p0 = theNeutron.GetTotalMomentum();
      G4double mN = theNeutron.GetMass();
      G4double mT = theTarget.GetMass();
      G4double eE = e0+mT;
      G4double ap = (mT+eE)*(mT-eE) + (p0+mN)*(p0-mN);
      G4double a = 4*(eE+p0*cosTh)*(eE-p0*cosTh);
      G4double b = 4*ap*p0*cosTh;
      G4double c = (2.*eE*mN-ap)*(2.*eE*mN+ap);
      G4double en = (-b+sqrt(b*b - 4*a*c) )/(2*a);
      G4ThreeVector tempVector(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
      theNeutron.SetMomentum(tempVector);
      theNeutron.SetTotalEnergy(sqrt(en*en+theNeutron.GetMass()*theNeutron.GetMass()));
      // first to lab     
      theNeutron.Lorentz(theNeutron, -1.*theTarget);
      // now to CMS
      theNeutron.Lorentz(theNeutron, theCMS);
      theTarget.SetMomentum(-theNeutron.GetMomentum());
      theTarget.SetTotalEnergy(theNeutron.GetTotalEnergy());
      // and back to lab
      theNeutron.Lorentz(theNeutron, -1.*theCMS);
      theTarget.Lorentz(theTarget, -1.*theCMS);      
    }
    else if (frameFlag == 2) // CMS
    {
      theNeutron.Lorentz(theNeutron, theCMS);
      theTarget.Lorentz(theTarget, theCMS);
      G4double en = theNeutron.GetTotalMomentum();
      G4ThreeVector tempVector(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
      theNeutron.SetMomentum(tempVector);
      theTarget.SetMomentum(-tempVector);
      G4double tP = theTarget.GetTotalMomentum();
      G4double tM = theTarget.GetMass();
      theTarget.SetTotalEnergy(sqrt((tP+tM)*(tP+tM)-2.*tP*tM));
      theNeutron.Lorentz(theNeutron, -1.*theCMS);
      theTarget.Lorentz(theTarget, -1.*theCMS);
    }
    else
    {
      G4cout <<"Value of frameFlag (1=LAB, 2=CMS): "<<frameFlag;
      G4Exception("G4NeutronHPElasticFS::ApplyYourSelf frameflag incorrect");
    }
    // now all in Lab
// nun den recoil generieren...und energy change, momentum change angeben.
    theResult.SetEnergyChange(theNeutron.GetKineticEnergy());
    theResult.SetMomentumChange(theNeutron.GetMomentum().unit());
    G4DynamicParticle* theRecoil = new G4DynamicParticle;
    if(targetMass<4.5)
    {
      G4bool He3flag = false;
      if(targetMass<1)
      {
        // proton
        theRecoil->SetDefinition(G4Proton::Proton());
      }
      else if(targetMass<2 )
      {
        // deuteron
        theRecoil->SetDefinition(G4Deuteron::Deuteron());
      }
      else if(targetMass<2.999 )
      {
        // 3He 
        theRecoil->SetDefinition(G4He3::He3());
      }
      else if(targetMass<3 )
      {
        // Triton
        theRecoil->SetDefinition(G4Triton::Triton());
      }
      else
      {
        // alpha
        theRecoil->SetDefinition(G4Alpha::Alpha());
      }
    }
    else
    {
      theRecoil->SetDefinition(G4ParticleTable::GetParticleTable()->FindIon(theBaseZ, theBaseA, 0, theBaseZ));
    }
    theRecoil->SetMomentum(theTarget.GetMomentum());
    theResult.SetNumberOfSecondaries(1);
    theResult.AddSecondary(theRecoil);
//    G4cout << "G4NeutronHPElasticFS::ApplyYourself 10+"<<G4endl;
    // postpone the tracking of the primary neutron
     theResult.SetStatusChange(fSuspend);
    return &theResult;
  }
