// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPEnAngCorrelation.hh"

  G4NeutronHPEnAngCorrelation::G4NeutronHPEnAngCorrelation()
  {
    theProducts = NULL;
    inCharge = false;
    theTotalMeanEnergy = -1.;
  }
  G4NeutronHPEnAngCorrelation::~G4NeutronHPEnAngCorrelation()
  {
    if(theProducts!=NULL) delete [] theProducts;
  }

G4ReactionProduct * G4NeutronHPEnAngCorrelation::SampleOne(G4double anEnergy)
{  
  G4ReactionProduct * result = new G4ReactionProduct;
  
  // do we have an appropriate distribution
  if(nProducts!=1) G4Exception("More than one product in SampleOne");
  
  // get the result
  G4ReactionProductVector * temp=NULL;
  G4int i=0;
  while(temp == NULL) temp = theProducts[i++].Sample(anEnergy);
  
  // is the multiplicity correct
  if(temp->length()!=1) G4Exception("SampleOne: Yield not correct");
  
  // fill result
  result = temp->at(0);
  
  // some garbage collection
  delete temp;
  
  // return result
  return result;
}

G4ReactionProductVector * G4NeutronHPEnAngCorrelation::Sample(G4double anEnergy)
{
  G4ReactionProductVector * result = new G4ReactionProductVector;
  G4int i;
  G4ReactionProductVector * it;
  G4ReactionProduct theCMS;
  if(frameFlag==2)
  {
    // simplify and double check @
    G4ThreeVector the3Neutron = theNeutron.GetMomentum();
    G4double nEnergy = theNeutron.GetTotalEnergy();
    G4ThreeVector the3Target = theTarget.GetMomentum();
    G4double tEnergy = theTarget.GetTotalEnergy();
    G4double totE = nEnergy+tEnergy;
    G4ThreeVector the3CMS = the3Target+the3Neutron;
    theCMS.SetMomentum(the3CMS);
    G4double cmsMom = sqrt(the3CMS*the3CMS);
    G4double sqrts = sqrt((totE-cmsMom)*(totE+cmsMom));
    theCMS.SetMass(sqrts);
    theCMS.SetTotalEnergy(totE);
    G4ReactionProduct aNeutron;
    aNeutron.Lorentz(theNeutron, theCMS);
    anEnergy = aNeutron.GetKineticEnergy();
  }
  theTotalMeanEnergy=0;
  for(i=0; i<nProducts; i++)
  {
    it = theProducts[i].Sample(anEnergy);
    G4double aMeanEnergy = theProducts[i].MeanEnergyOfThisInteraction();
    if(aMeanEnergy>0)
    {
      theTotalMeanEnergy += aMeanEnergy;
    }
    else
    {
      theTotalMeanEnergy = anEnergy/nProducts+theProducts[i].GetQValue();
    }
    if(it!=NULL)
    {
      for(G4int ii=0; ii<it->length(); ii++)
      {
	if(frameFlag==1) // target rest
	{
          it->at(ii)->Lorentz(*(it->at(ii)), -1.*theTarget);
	}
	else if(frameFlag==2) // CMS
	{
#ifdef G4_NHP_DEBUG
          cout <<"G4NeutronHPEnAngCorrelation: "<<
        	 it->at(ii)->GetTotalEnergy()<<" "<<
        	 it->at(ii)->GetMomentum()<<endl;
#endif
          it->at(ii)->Lorentz(*(it->at(ii)), -1.*theCMS);
	}
	else
	{
          G4Exception("G4NeutronHPEnAngCorrelation::Sample: The frame of the finalstate is not specified");
	}
	result->insert(it->at(ii));
      }
    delete it;
    }
  }   
  return result;
}
