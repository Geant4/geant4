// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPAngular.hh"

void G4NeutronHPAngular::Init(ifstream & aDataFile)
{
//  G4cout << "here we are entering the Angular Init"<<endl;
  G4int dummy;
  aDataFile >> theAngularDistributionType >> targetMass;
  aDataFile >> frameFlag;
  if(theAngularDistributionType == 0)
  {
    theIsoFlag = true;
  }
  else if(theAngularDistributionType==1)
  {
    G4int nEnergy;
    aDataFile >> nEnergy;  
    theCoefficients = new G4NeutronHPLegendreStore(nEnergy);
    theCoefficients->InitInterpolation(aDataFile);
    G4double temp, energy;
    G4int tempdep, nLegendre;
    G4int i, ii;
    for (i=0; i<nEnergy; i++)
    {
      aDataFile >> temp >> energy >> tempdep >> nLegendre;
      energy *=eV;
      theCoefficients->Init(i, energy, nLegendre);
      theCoefficients->SetTemperature(i, temp);
      G4double coeff=0;
      for(ii=0; ii<nLegendre; ii++)
      {
        aDataFile >> coeff;
        theCoefficients->SetCoeff(i, ii+1, coeff);
      }
    }
  }
  else if (theAngularDistributionType==2)
  {
    G4int nEnergy;
    aDataFile >> nEnergy;
    theProbArray = new G4NeutronHPPartial(nEnergy, nEnergy);
    theProbArray->InitInterpolation(aDataFile);
    G4double temp, energy;
    G4int tempdep, nPoints;
    for(G4int i=0; i<nEnergy; i++)
    {
      aDataFile >> temp >> energy >> tempdep;
      energy *= eV;
      theProbArray->SetT(i, temp);
      theProbArray->SetX(i, energy);
      theProbArray->InitData(i, aDataFile);
    }
  }
  else
  {
    theIsoFlag = false;
    G4cout << "unknown distribution found for Angular"<<endl;
    G4Exception("unknown distribution needs implementation!!!");
  }    
}

void G4NeutronHPAngular::SampleAndUpdate(G4ReactionProduct & aHadron)
{
  if(theIsoFlag)
  {
//  G4cout << "Angular result "<<aHadron.GetTotalMomentum()<<" ";
      G4double costheta = 2.*G4UniformRand()-1;
      G4double theta = acos(costheta);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = sin(theta);
      G4double en = aHadron.GetTotalMomentum();
      G4ThreeVector temp(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
      aHadron.SetMomentum( temp );
      aHadron.Lorentz(aHadron, -1.*theTarget);
  }
  else
  {
    if(theAngularDistributionType == 1) // LAB
    {
      G4double en = aHadron.GetTotalMomentum();
      G4ReactionProduct boosted;
      boosted.Lorentz(theNeutron, theTarget);
      G4double kineticEnergy = boosted.GetKineticEnergy();
      G4double cosTh = theCoefficients->SampleMax(kineticEnergy);
      G4double theta = acos(cosTh);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = sin(theta);
      G4ThreeVector temp(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
      aHadron.SetMomentum( temp );
    }
    else if(theAngularDistributionType == 2) // costh in CMS
    {
      G4ReactionProduct boostedN;
      boostedN.Lorentz(theNeutron, theTarget);
      G4double kineticEnergy = boostedN.GetKineticEnergy();
      G4double cosTh = theProbArray->Sample(kineticEnergy); 
      G4double theta = acos(cosTh);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = sin(theta);      
      
      G4ThreeVector temp(sinth*cos(phi), sinth*sin(phi), cos(theta) ); //CMS
      G4double en = aHadron.GetTotalEnergy(); // Target rest
      
      // get trafo from Target rest frame to CMS
      G4ReactionProduct boostedT;
      boostedT.Lorentz(theTarget, theTarget);
      
      G4ThreeVector the3Neutron = boostedN.GetMomentum();
      G4double nEnergy = boostedN.GetTotalEnergy();
      G4ThreeVector the3Target = boostedT.GetMomentum();
      G4double tEnergy = boostedT.GetTotalEnergy();
      G4double totE = nEnergy+tEnergy;
      G4ThreeVector the3trafo = -the3Target-the3Neutron;
      G4ReactionProduct trafo; // for transformation from CMS to target rest frame
      trafo.SetMomentum(the3trafo);
      G4double cmsMom = sqrt(the3trafo*the3trafo);
      G4double sqrts = sqrt((totE-cmsMom)*(totE+cmsMom));
      trafo.SetMass(sqrts);
      trafo.SetTotalEnergy(totE);
      
      G4double gamma = trafo.GetTotalEnergy()/trafo.GetMass();
      G4double cosalpha = temp*trafo.GetMomentum()/trafo.GetTotalMomentum()/temp.mag();
      G4double fac = cosalpha*trafo.GetTotalMomentum()/trafo.GetMass();
      fac*=gamma;
      
      G4double mom;
      mom = sqrt( en*fac*en*fac - 
                   (fac*fac - gamma*gamma)*
                   (en*en - gamma*gamma*aHadron.GetMass()*aHadron.GetMass())
                );
      mom = -en*fac - mom;
      mom /= (fac*fac-gamma*gamma);
      temp = mom*temp;
      
      aHadron.SetMomentum( temp ); // now all in CMS
      aHadron.SetTotalEnergy( sqrt( mom*mom + aHadron.GetMass()*aHadron.GetMass() ) );
      aHadron.Lorentz(aHadron, trafo); // now in target rest frame
    }
    else
    {
      G4Exception("Tried to sample non isotropic neutron angular");
    }
  }
  aHadron.Lorentz(aHadron, -1.*theTarget); 
//  G4cout << aHadron.GetMomentum()<<" ";
//  G4cout << aHadron.GetTotalMomentum()<<endl;
}
