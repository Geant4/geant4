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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 080801 Prohibit level transition to oneself by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPDeExGammas.hh"
#include "G4SystemOfUnits.hh"

void G4ParticleHPDeExGammas::Init(std::istream & aDataFile)
{
  //  G4cout << this << "ExGammas Init LEVEL " << G4endl; //GDEB
  G4ParticleHPGamma ** theGammas = new G4ParticleHPGamma * [50];
  G4int nGammas = 0;
  G4int nBuff = 50;
  for(;;)
  {
    G4ParticleHPGamma * theNew = new G4ParticleHPGamma;
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
        G4ParticleHPGamma ** buffer = new G4ParticleHPGamma * [nBuff];
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
    if(std::abs(currentE-nextE)>epsilon) nLevels++;
    currentE = nextE;
  }

  //  G4cout << this << "LEVEL " << nLevels << G4endl; //GDEB
  // Build the levels

  theLevels = new G4ParticleHPLevel[nLevels];
  levelStart = new G4int [nLevels];
  levelSize = new G4int [nLevels];

  // fill the levels

  currentE = 0;
  nextE = 0;
  G4int levelCounter=-1;
  for(i=0; i<nGammas; i++)
  {
    nextE = theGammas[i]->GetLevelEnergy();
    if(std::abs(currentE-nextE)>epsilon) 
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
      if(std::abs(currentLevelE-(levelE+gammaE))<min)
      {
        min = std::abs(currentLevelE-(levelE+gammaE));
        it = ii;
      }
    }
//080728
    if ( it != -1 && currentLevelE == theLevels[it].GetLevelEnergy() )
    {
       //TK Comment; Some data file in /Inelastic/Gammas has inconsistent level data (no level to transit)
       //G4cout << "DeExGammas Transition level error: it " << it << " " << currentLevelE << " " << gammaE << " " << theLevels[it-1].GetLevelEnergy() << " " << currentLevelE - theLevels[it-1].GetLevelEnergy() << G4endl;
       // Forced to connect the next(previous) level 
       it +=-1;
    }
//080728
    if(it!=-1) theGammas[i]->SetNext(&theLevels[it]);
  }
  // some garbage collection

  delete [] theGammas;

  // and we are Done.
}
