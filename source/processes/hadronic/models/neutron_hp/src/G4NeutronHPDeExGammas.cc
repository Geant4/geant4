// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPDeExGammas.hh"

void G4NeutronHPDeExGammas::Init(ifstream & aDataFile)
{
  G4NeutronHPGamma ** theGammas = new G4NeutronHPGamma * [50];
  G4int nGammas = 0;
  G4int nBuff = 50;
  for(;;)
  {
    G4NeutronHPGamma * theNew = new G4NeutronHPGamma;
    if(!theNew->Init(aDataFile))
    {
      delete theNew;
      break;
    }
    else
    {
      if(nGammas==nBuff)
      {
        nBuff+=50;
        G4NeutronHPGamma ** buffer = new G4NeutronHPGamma * [nBuff];
        for(G4int i=0;i<nGammas;i++) buffer[i] = theGammas[i];
        delete [] theGammas;
        theGammas = buffer;
      }
      theGammas[nGammas] = theNew;
      nGammas++;
    }
  }
  // all gammas are in. Now sort them into levels.

  // count the levels

  G4double currentE = 0;
  G4double nextE = 0;
  G4int i;
  G4double epsilon = 0.01*keV;
  for(i=0; i<nGammas; i++)
  {
    nextE = theGammas[i]->GetLevelEnergy();
    if(abs(currentE-nextE)>epsilon) nLevels++;
    currentE = nextE;
  }

  // Build the levels

  theLevels = new G4NeutronHPLevel[nLevels];
  levelStart = new G4int [nLevels];
  levelSize = new G4int [nLevels];

  // fill the levels

  currentE = 0;
  nextE = 0;
  G4int levelCounter=-1;
  for(i=0; i<nGammas; i++)
  {
    nextE = theGammas[i]->GetLevelEnergy();
    if(abs(currentE-nextE)>epsilon) 
    {
      levelCounter++;
      levelStart[levelCounter] = i;
      levelSize[levelCounter] = 0;
    }
    levelSize[levelCounter]++;
    currentE = nextE;
  }

  for(i=0; i<nLevels; i++)
  {
    theLevels[i].SetNumberOfGammas(levelSize[i]);
    for(G4int ii=levelStart[i]; ii<levelStart[i]+levelSize[i]; ii++)
    {
      theLevels[i].SetGamma(ii-levelStart[i], theGammas[ii]);
    }
  }

// set the next relation in the gammas.
  G4double levelE, gammaE, currentLevelE;
  G4double min; 
  for(i=0; i<nGammas; i++)
  {
    G4int it=-1;
    gammaE = theGammas[i]->GetGammaEnergy();
    currentLevelE = theGammas[i]->GetLevelEnergy();
    min = currentLevelE-gammaE-epsilon;
    for(G4int ii=0; ii<nLevels; ii++)
    {
      levelE = theLevels[ii].GetLevelEnergy();
      if(abs(currentLevelE-(levelE+gammaE))<min)
      {
        min = abs(currentLevelE-(levelE+gammaE));
        it = ii;
      }
    }
    if(it!=-1) theGammas[i]->SetNext(&theLevels[it]);
  }
  // some garbage collection

  delete [] theGammas;

  // and we are Done.
}
