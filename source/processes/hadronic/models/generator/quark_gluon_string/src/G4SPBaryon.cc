#include "G4SPBaryon.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"

void G4SPBaryon::
SampleQuarkAndDiquark(G4int & quark, G4int & diQuark) const
{
  G4double random = 0.999999*G4UniformRand();
  G4double sum = 0;
  for(G4int i=0; i<thePartonInfo.length(); i++)
  {
    sum += thePartonInfo[i]->GetProbability();
    if (sum > random) 
    {
      if (theDefinition->GetPDGEncoding() < 0)
      {
	quark = -thePartonInfo[i]->GetDiQuark();
	diQuark = -thePartonInfo[i]->GetQuark();                   
      }
      else
      {
	quark = thePartonInfo[i]->GetQuark();
	diQuark = thePartonInfo[i]->GetDiQuark();                   
      }
      break;   
    }    
  }
}

void G4SPBaryon::
FindDiquark(G4int quark, G4int & diQuark) const
{
  G4double sum = 0;
  G4int i;
  for(i = 0; i<thePartonInfo.length(); i++)
  {
    if (abs(thePartonInfo[i]->GetQuark()) == abs(quark))
    {
      sum += thePartonInfo[i]->GetProbability();
    }
  }
  G4double random = 0.999999*G4UniformRand();
  G4double running = 0;
  for(i = 0; thePartonInfo[i]->GetQuark() != 0 && i < 6; i++)
  {
    if (abs(thePartonInfo[i]->GetQuark()) == abs(quark))
    {
      running += thePartonInfo[i]->GetProbability();
      if (running/sum >= random) 
      {
        diQuark = thePartonInfo[i]->GetDiQuark(); 
        break;                  
      }
    }
  }
}


G4SPBaryon::
G4SPBaryon(G4Proton * aProton)
{
theDefinition = aProton;
thePartonInfo.insert(new G4SPPartonInfo(2203, 1, 1./3.));
thePartonInfo.insert(new G4SPPartonInfo(2103, 2, 1./6.));
thePartonInfo.insert(new G4SPPartonInfo(2101, 2, 1./2.));
}

G4SPBaryon::
G4SPBaryon(G4Neutron * aNeutron)
{
theDefinition = aNeutron;
thePartonInfo.insert(new G4SPPartonInfo(2103, 1, 1./6.));
thePartonInfo.insert(new G4SPPartonInfo(2101, 1, 1./2.));
thePartonInfo.insert(new G4SPPartonInfo(1103, 2, 1./3.));
}

G4SPBaryon::
G4SPBaryon(G4SigmaPlus * aSigmaPlus)
{
theDefinition = aSigmaPlus;
thePartonInfo.insert(new G4SPPartonInfo(2203, 3, 1./3.));
thePartonInfo.insert(new G4SPPartonInfo(3203, 2, 1./6.));
thePartonInfo.insert(new G4SPPartonInfo(3201, 2, 1./2.));
}

G4SPBaryon::
G4SPBaryon(G4SigmaZero * aSigmaZero)
{
theDefinition = aSigmaZero;
thePartonInfo.insert(new G4SPPartonInfo(2103, 3, 1./3.));
thePartonInfo.insert(new G4SPPartonInfo(3203, 1, 1./12.));
thePartonInfo.insert(new G4SPPartonInfo(3201, 1, 1./4.));
thePartonInfo.insert(new G4SPPartonInfo(3103, 2, 1./12.));
thePartonInfo.insert(new G4SPPartonInfo(3101, 2, 1./4.));
}

G4SPBaryon::
G4SPBaryon(G4SigmaMinus * aSigmaMinus)
{
theDefinition = aSigmaMinus;
thePartonInfo.insert(new G4SPPartonInfo(1103, 3, 1./3.));
thePartonInfo.insert(new G4SPPartonInfo(3103, 1, 1./6.));
thePartonInfo.insert(new G4SPPartonInfo(3101, 1, 1./2.));
}

G4SPBaryon::
G4SPBaryon(G4XiMinus * aXiMinus)
{
theDefinition = aXiMinus;
thePartonInfo.insert(new G4SPPartonInfo(3103, 3, 1./6.));
thePartonInfo.insert(new G4SPPartonInfo(3101, 3, 1./2.));
thePartonInfo.insert(new G4SPPartonInfo(3303, 1, 1./3.));
}
G4SPBaryon::
G4SPBaryon(G4XiZero * aXiZero)
{
theDefinition = aXiZero;
thePartonInfo.insert(new G4SPPartonInfo(3203, 3, 1./6.));
thePartonInfo.insert(new G4SPPartonInfo(3201, 3, 1./2.));
thePartonInfo.insert(new G4SPPartonInfo(3303, 2, 1./3.));
}
G4SPBaryon::
G4SPBaryon(G4Lambda * aLambda)
{
theDefinition = aLambda;
thePartonInfo.insert(new G4SPPartonInfo(2103, 3, 1./3.));
thePartonInfo.insert(new G4SPPartonInfo(3203, 1, 1./4.));
thePartonInfo.insert(new G4SPPartonInfo(3201, 1, 1./12.));
thePartonInfo.insert(new G4SPPartonInfo(3103, 2, 1./4.));
thePartonInfo.insert(new G4SPPartonInfo(3101, 2, 1./12.));
}
G4SPBaryon::
G4SPBaryon(G4OmegaMinus * anOmegaMinus)
{
theDefinition = anOmegaMinus;
thePartonInfo.insert(new G4SPPartonInfo(3303, 3, 1.));
}
// non static particles
G4SPBaryon::
G4SPBaryon(G4ParticleDefinition * aDefinition)
{
  theDefinition = aDefinition;
  if(theDefinition == 
     G4ParticleTable::GetParticleTable()->FindParticle(2224))// D++
  {
    thePartonInfo.insert(new G4SPPartonInfo(2203, 2, 1.));
  }
  else if(theDefinition == 
     G4ParticleTable::GetParticleTable()->FindParticle(2214))// D+
  {
    thePartonInfo.insert(new G4SPPartonInfo(2203, 1, 1./3.));
    thePartonInfo.insert(new G4SPPartonInfo(2103, 2, 2./3.));
  }
  else if(theDefinition ==
     G4ParticleTable::GetParticleTable()->FindParticle(2114))// D0
  {
    thePartonInfo.insert(new G4SPPartonInfo(2103, 1, 2./3.));
    thePartonInfo.insert(new G4SPPartonInfo(2103, 2, 1./3.));
  }
  else if(theDefinition == 
     G4ParticleTable::GetParticleTable()->FindParticle(1114))// D-
  {
    thePartonInfo.insert(new G4SPPartonInfo(1103, 1, 1.));
  }
  else if(theDefinition == 
     G4ParticleTable::GetParticleTable()->FindParticle(3224))// S*+
  {
    thePartonInfo.insert(new G4SPPartonInfo(2203, 3, 1./3.));
    thePartonInfo.insert(new G4SPPartonInfo(3203, 2, 2./3.));
  }
  else if(theDefinition == 
     G4ParticleTable::GetParticleTable()->FindParticle(3214))// S*0
  {
    thePartonInfo.insert(new G4SPPartonInfo(2103, 3, 1./3.));
    thePartonInfo.insert(new G4SPPartonInfo(3203, 1, 1./3.));
    thePartonInfo.insert(new G4SPPartonInfo(3103, 2, 1./3.));
  }
  else if(theDefinition == 
     G4ParticleTable::GetParticleTable()->FindParticle(3224))// S*-
  {
    thePartonInfo.insert(new G4SPPartonInfo(1103, 3, 1./3.));
    thePartonInfo.insert(new G4SPPartonInfo(3103, 1, 2./3.));
  }
  else if(theDefinition == 
     G4ParticleTable::GetParticleTable()->FindParticle(3324))// Xi*0
  {
    thePartonInfo.insert(new G4SPPartonInfo(3203, 3, 1./3.));
    thePartonInfo.insert(new G4SPPartonInfo(3303, 2, 2./3.));
  }
  else if(theDefinition == 
     G4ParticleTable::GetParticleTable()->FindParticle(3314))// Xi*-
  {
    thePartonInfo.insert(new G4SPPartonInfo(3103, 3, 2./3.));
    thePartonInfo.insert(new G4SPPartonInfo(3303, 1, 1./3.));
  }
}
