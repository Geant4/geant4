// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// there is a lot of unused (and undebugged) code in this file. Kept for the moment just in case. @@

#include "G4NeutronHPPhotonDist.hh"
#include "G4NeutronHPLegendreStore.hh"
#include "G4Electron.hh"

  G4NeutronHPPhotonDist::G4NeutronHPPhotonDist()
  {
     disType = NULL;
     energy = NULL;
     theYield = NULL;
     thePartialXsec = NULL;
     isPrimary = NULL;
     theShells = NULL;
     theGammas = NULL;
     nNeu = NULL;
     theLegendre = NULL;
     theAngular = NULL;
     distribution = NULL;
     probs = NULL;
     partials = NULL;
     actualMult = NULL;

     theLevelEnergies = NULL;
     theTransitionProbabilities = NULL;
     thePhotonTransitionFraction = NULL;
  }

  G4NeutronHPPhotonDist::~G4NeutronHPPhotonDist()
  {
     if(disType != NULL) delete [] disType;
     if(energy != NULL) delete [] energy;
     if(theYield != NULL) delete [] theYield;
     if(thePartialXsec != NULL) delete [] thePartialXsec;
     if(isPrimary != NULL) delete [] isPrimary;
     if(theShells != NULL) delete [] theShells;
     if(theGammas != NULL) delete [] theGammas;
     if(nNeu != NULL) delete [] nNeu;
     if(theLegendre != NULL) delete [] theLegendre;
     if(theAngular != NULL) delete [] theAngular;
     if(distribution != NULL) delete [] distribution;
     if(probs != NULL) delete [] probs;
     if(partials != NULL) delete [] partials;
     if(actualMult != NULL) delete [] actualMult;

     if(theLevelEnergies != NULL) delete theLevelEnergies;
     if(theTransitionProbabilities != NULL) delete theTransitionProbabilities;
     if(thePhotonTransitionFraction != NULL) delete thePhotonTransitionFraction;
  }

G4bool G4NeutronHPPhotonDist::InitMean(ifstream & aDataFile)
{
  G4bool result = true;
  if(aDataFile >> repFlag)
  {
    aDataFile >> targetMass;
    G4int d1=0;
    G4double e=0, y=0, z=0;
    if(repFlag==1)
    {
    // multiplicities
      aDataFile >> nDiscrete;
      disType = new G4int[nDiscrete];
      energy = new G4double[nDiscrete];
      actualMult = new G4int[nDiscrete];
      theYield = new G4NeutronHPVector[nDiscrete];
      for (G4int i=0; i<nDiscrete; i++)
      {
	aDataFile >> disType[i]>>energy[i];
	energy[i]*=eV;
	theYield[i].Init(aDataFile, eV);
      }
    }
    else if(repFlag == 2)
    {
       aDataFile >> theInternalConversionFlag;
       aDataFile >> theBaseEnergy;
       theBaseEnergy*=eV;
       aDataFile >> theInternalConversionFlag;
       aDataFile >> nGammaEnergies;
       theLevelEnergies = new G4double[nGammaEnergies];
       theTransitionProbabilities = new G4double[nGammaEnergies];
       if(theInternalConversionFlag == 2) thePhotonTransitionFraction = new G4double[nGammaEnergies];
       for(G4int  ii=0; ii<nGammaEnergies; ii++)
       {
	 if(theInternalConversionFlag == 1)
	 {
           aDataFile >> theLevelEnergies[ii] >> theTransitionProbabilities[ii];
	   theLevelEnergies[ii]*=eV;
	 }
	 else if(theInternalConversionFlag == 2)
	 {
           aDataFile >> theLevelEnergies[ii] >> theTransitionProbabilities[ii] >> thePhotonTransitionFraction[ii];
	   theLevelEnergies[ii]*=eV;
	 }
	 else
	 {
           G4Exception("G4NeutronHPPhotonDist: Unknown conversion flag");
	 }
      }
       // Note, that this is equivalent to using the 'Gamma' classes.
      // G4Exception("G4NeutronHPPhotonDist: Transition probability array not sampled for the moment.");
    }
    else
    {
      G4cout << "Data representation in G4NeutronHPPhotonDist: "<<repFlag<<endl;
      G4Exception("G4NeutronHPPhotonDist: This data representation is not implemented.");
    }
  }
  else
  {
    result = false;
  }
  return result;
}

void G4NeutronHPPhotonDist::InitAngular(ifstream & aDataFile)
{
  G4int i, ii;
  //angular distributions
  aDataFile >> isoFlag;
  if (isoFlag != 1)
  {
    aDataFile >> tabulationType >> nDiscrete2 >> nIso;
    theShells = new G4double[nDiscrete2];
    theGammas = new G4double[nDiscrete2];
    for (i=0; i< nIso; i++) // isotropic photons
    {
        aDataFile >> theGammas[i] >> theShells[i];
        theGammas[i]*=eV;
        theShells[i]*=eV;
    }
    G4double eNeu, coeff;
    G4int nPoly, nProb;
    nNeu = new G4int [nDiscrete2-nIso];
    if(tabulationType==1)theLegendre=new G4NeutronHPLegendreTable *[nDiscrete2-nIso];
    if(tabulationType==2)theAngular =new G4NeutronHPAngularP *[nDiscrete2-nIso];
    for(i=nIso; i< nDiscrete2; i++)
    {
      if(tabulationType==1) 
      {
        aDataFile >> theGammas[i] >> theShells[i] >> nNeu[i-nIso];
        theGammas[i]*=eV;
        theShells[i]*=eV;
        theLegendre[i-nIso]=new G4NeutronHPLegendreTable[nNeu[i-nIso]];
        theLegendreManager.Init(aDataFile); 
        for (ii=0; ii<nNeu[i-nIso]; ii++)
        {
          theLegendre[i-nIso][ii].Init(aDataFile);
        }
      }
      else if(tabulationType==2)
      {
        aDataFile >> theGammas[i] >> theShells[i] >> nNeu[i-nIso];
        theGammas[i]*=eV;
        theShells[i]*=eV;
        theAngular[i-nIso]=new G4NeutronHPAngularP[nNeu[i-nIso]];
        for (ii=0; ii<nNeu[i-nIso]; ii++)
        {
          theAngular[i-nIso][ii].Init(aDataFile);
        }
      }
      else
      {
        G4cout << "tabulation type: tabulationType"<<endl;
        G4Exception("cannot deal with this tabulation type for angular distributions.");
      }
    }
  }
}


void G4NeutronHPPhotonDist::InitEnergies(ifstream & aDataFile)
{
  G4int i, energyDistributionsNeeded = 0;
  for (i=0; i<nDiscrete; i++)
  {
    if( disType[i]==1) energyDistributionsNeeded =1;
  }
  if(!energyDistributionsNeeded) return;
  aDataFile >>  nPartials;
  distribution = new G4int[nPartials];
  probs = new G4NeutronHPVector[nPartials];
  partials = new G4NeutronHPPartial * [nPartials];
  G4int nen;
  G4int dummy;
  for (i=0; i<nPartials; i++)
  {
    aDataFile >> dummy;
    probs[i].Init(aDataFile, eV);
    aDataFile >> nen;
    partials[i] = new G4NeutronHPPartial(nen);
    partials[i]->InitInterpolation(aDataFile);
    partials[i]->Init(aDataFile);
  }
}

void G4NeutronHPPhotonDist::InitPartials(ifstream & aDataFile)
{
  aDataFile >> nDiscrete >> targetMass;
  if(nDiscrete != 1)
  {
    theTotalXsec.Init(aDataFile, eV);
  }
  G4int i;
  theGammas = new G4double[nDiscrete];
  theShells = new G4double[nDiscrete];
  isPrimary = new G4int[nDiscrete];
  disType = new G4int[nDiscrete];
  thePartialXsec = new G4NeutronHPVector[nDiscrete];
  for(i=0; i<nDiscrete; i++)
  {
    aDataFile>>theGammas[i]>>theShells[i]>>isPrimary[i]>>disType[i];
    theGammas[i]*=eV;
    theShells[i]*=eV;
    thePartialXsec[i].Init(aDataFile, eV);
  }  
}

G4ReactionProductVector * G4NeutronHPPhotonDist::GetPhotons(G4double anEnergy)
{
  // the partial cross-section case is not in this yet. @@@@
  G4int i, ii, iii;
  G4int nSecondaries = 0;
  G4ReactionProductVector * thePhotons = new G4ReactionProductVector;
  if(repFlag==1)
  {
    G4double current=0;
    for(i=0; i<nDiscrete; i++)
    {
      current = theYield[i].GetY(anEnergy);
      actualMult[i] = RandPoisson::shoot(current); // max cut-off still missing @@@
      if(nDiscrete==1&&current<1.0001) 
      {
        actualMult[i] = current;
        if(current<1) 
        {
          actualMult[i] = 0;
          if(G4UniformRand()<current) actualMult[i] = 1;
        }
      }
      nSecondaries += actualMult[i];
    }
    for(i=0;i<nSecondaries;i++)
    {
      G4ReactionProduct * theOne = new G4ReactionProduct;
      theOne->SetDefinition(G4Gamma::Gamma());
      thePhotons->insert(theOne);
    }
    G4int count=0;
    for(i=0; i<nDiscrete; i++)
    { 
      for(ii=0; ii< actualMult[i]; ii++)
      {   
	if(disType[i]==1) // continuum
	{
          G4double econt=0, sum=0, run=0;
          for(iii=0; iii<nPartials; iii++) sum+=probs[iii].GetY(anEnergy);
          G4double random = G4UniformRand();
          G4int theP = 0;
          for(iii=0; iii<nPartials; iii++)
          {
            run+=probs[iii].GetY(anEnergy);
            theP = iii;
            if(random<run/sum) break;
          }
          if(theP==nPartials) theP=nPartials-1; // das sortiert J aus.
          sum=0; 
          G4NeutronHPVector * temp;
          temp = partials[theP]->GetY(anEnergy); //@@@ look at, seems fishy
          G4double eGamm = temp->Sample();
          thePhotons->at(count)->SetKineticEnergy(eGamm);
          delete temp;
	}
	else // discrete
	{
          thePhotons->at(count)->SetKineticEnergy(energy[i]);
	}
	count++;
	if(count > nSecondaries)  G4Exception("G4NeutronHPPhotonDist::GetPhotons inconsistancy");
      }
    }
    // now do the angular distributions...
    G4double count1=0;
    if( isoFlag == 1)
    {
      for (i=0; i< nSecondaries; i++)
      {
	G4double costheta = 2.*G4UniformRand()-1;
	G4double theta = acos(costheta);
	G4double phi = twopi*G4UniformRand();
	G4double sinth = sin(theta);
	G4double en = thePhotons->at(i)->GetTotalEnergy();
	G4ThreeVector temp(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
	thePhotons->at(i)->SetMomentum( temp ) ;
  //      G4cout << "Isotropic distribution in PhotonDist"<<temp<<endl;
      }
    }
    else
    {
      for(i=0; i<nSecondaries; i++)
      { 
	G4double currentEnergy = thePhotons->at(i)->GetTotalEnergy();
	for(ii=0; ii<nDiscrete2; ii++) 
	{
          if (abs(currentEnergy-theGammas[ii])<0.1*keV) break;
	}
	if(ii==nDiscrete2) ii--; // fix for what seems an (file12 vs file 14) inconsistancy found in the ENDF 7N14 data. @@
	if(ii<nIso)
	{
          // isotropic distribution
          G4double theta = pi*G4UniformRand();
          G4double phi = twopi*G4UniformRand();
          G4double sinth = sin(theta);
          G4double en = thePhotons->at(i)->GetTotalEnergy();
          G4ThreeVector tempVector(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
          thePhotons->at(i)->SetMomentum( tempVector ) ;
	}
	else if(tabulationType==1)
	{
          // legendre polynomials
          G4double rand = G4UniformRand();
          G4int it;
          for (iii=0; iii<nNeu[ii-nIso]; iii++) // find the neutron energy
          {
            it = iii;
	    if(theLegendre[ii-nIso][iii].GetEnergy()>anEnergy)
              break;
          }
          G4NeutronHPLegendreStore aStore(2);
          aStore.SetCoeff(1, &(theLegendre[ii-nIso][it]));  
          aStore.SetCoeff(0, &(theLegendre[ii-nIso][it-1])); 
          G4double cosTh = aStore.SampleMax(anEnergy);
          G4double theta = acos(cosTh);
          G4double phi = twopi*G4UniformRand();
          G4double sinth = sin(theta);
          G4double en = thePhotons->at(i)->GetTotalEnergy();
          G4ThreeVector tempVector(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
          thePhotons->at(i)->SetMomentum( tempVector ) ;
	}
	else
	{
          // tabulation of probabilities.
          G4int it;
          for (iii=0; iii<nNeu[ii-nIso]; iii++) // find the neutron energy
          {
            it = iii;
	    if(theAngular[ii-nIso][iii].GetEnergy()>anEnergy)
              break;
          }
          G4double costh = theAngular[ii-nIso][it].GetCosTh(); // no interpolation yet @@
          G4double theta = acos(costh);
          G4double phi = twopi*G4UniformRand();
          G4double sinth = sin(theta);
          G4double en = thePhotons->at(i)->GetTotalEnergy();
          G4ThreeVector tmpVector(en*sinth*cos(phi), en*sinth*sin(phi), en*costh );
          thePhotons->at(i)->SetMomentum( tmpVector ) ;
	}
      }  
    } 
  }
  else if(repFlag == 2)
  {
    G4double * running = new G4double[nGammaEnergies];
    running[0]=theTransitionProbabilities[0];
    G4int i;
    for(i=1; i<nGammaEnergies; i++)
    {
      running[i]+=theTransitionProbabilities[i];
    }
    G4double random = G4UniformRand();
    G4int it=0;
    for(i=0; i<nGammaEnergies; i++)
    {
      it = i;
      if(random < running[i]/running[nGammaEnergies-1]) break;
    }
    delete [] running;
    G4double totalEnergy = theBaseEnergy - theLevelEnergies[it];
    G4ReactionProduct * theOne = new G4ReactionProduct;
    theOne->SetDefinition(G4Gamma::Gamma());
    random = G4UniformRand();
    if(theInternalConversionFlag==2 && random>thePhotonTransitionFraction[it])
    {
      theOne->SetDefinition(G4Electron::Electron());
    }
    theOne->SetTotalEnergy(totalEnergy);
    G4double count1=0;
    if( isoFlag == 1)
    {
      G4double costheta = 2.*G4UniformRand()-1;
      G4double theta = acos(costheta);
      G4double phi = twopi*G4UniformRand();
      G4double sinth = sin(theta);
      G4double en = theOne->GetTotalEnergy();
      G4ThreeVector temp(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
      theOne->SetMomentum( temp ) ;
    }
    else
    {
      G4double currentEnergy = theOne->GetTotalEnergy();
      for(ii=0; ii<nDiscrete2; ii++) 
      {
        if (abs(currentEnergy-theGammas[ii])<0.1*keV) break;
      }
      if(ii==nDiscrete2) ii--; // fix for what seems an (file12 vs file 14) inconsistancy found in the ENDF 7N14 data. @@
      if(ii<nIso)
      {
        // isotropic distribution
        G4double theta = pi*G4UniformRand();
        G4double phi = twopi*G4UniformRand();
        G4double sinth = sin(theta);
        G4double en = theOne->GetTotalEnergy();
        G4ThreeVector tempVector(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
        theOne->SetMomentum( tempVector ) ;
      }
      else if(tabulationType==1)
      {
        // legendre polynomials
        G4double rand = G4UniformRand();
        G4int it;
        for (iii=0; iii<nNeu[ii-nIso]; iii++) // find the neutron energy
        {
          it = iii;
	  if(theLegendre[ii-nIso][iii].GetEnergy()>anEnergy)
            break;
        }
        G4NeutronHPLegendreStore aStore(2);
        aStore.SetCoeff(1, &(theLegendre[ii-nIso][it]));  
        aStore.SetCoeff(0, &(theLegendre[ii-nIso][it-1])); 
        G4double cosTh = aStore.SampleMax(anEnergy);
        G4double theta = acos(cosTh);
        G4double phi = twopi*G4UniformRand();
        G4double sinth = sin(theta);
        G4double en = theOne->GetTotalEnergy();
        G4ThreeVector tempVector(en*sinth*cos(phi), en*sinth*sin(phi), en*cos(theta) );
        theOne->SetMomentum( tempVector ) ;
      }
      else
      {
        // tabulation of probabilities.
        G4int it;
        for (iii=0; iii<nNeu[ii-nIso]; iii++) // find the neutron energy
        {
          it = iii;
	  if(theAngular[ii-nIso][iii].GetEnergy()>anEnergy)
            break;
        }
        G4double costh = theAngular[ii-nIso][it].GetCosTh(); // no interpolation yet @@
        G4double theta = acos(costh);
        G4double phi = twopi*G4UniformRand();
        G4double sinth = sin(theta);
        G4double en = theOne->GetTotalEnergy();
        G4ThreeVector tmpVector(en*sinth*cos(phi), en*sinth*sin(phi), en*costh );
        theOne->SetMomentum( tmpVector ) ;
      }
    }
    thePhotons->insert(theOne);
  }
  else
  {
    delete thePhotons;
    thePhotons = NULL; // no gamma data available; some work needed @@@@@@@
  }    
  return thePhotons;
}

