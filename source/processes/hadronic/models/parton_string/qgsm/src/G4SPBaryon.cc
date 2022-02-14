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
#include "G4SPBaryon.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"

// correcting numbers, HPW Dec 1999

G4int G4SPBaryon::FindQuark(G4int diQuark) const
{
  G4double sum = GetProbability(diQuark);
  G4double random = G4UniformRand();
  G4double running = 0;
  G4int Quark(0);
  typedef std::vector<G4SPPartonInfo *>::const_iterator iter;
  iter i;
  for (i = thePartonInfo.begin(); i!=thePartonInfo.end(); i++)
  {
    if (std::abs((*i)->GetDiQuark()) == std::abs(diQuark))
    {
      running += (*i)->GetProbability();
      if (running/sum >= random)
      {
        Quark = (*i)->GetQuark();
	break;
      }
    }
  }
  return Quark;
}


G4double G4SPBaryon::GetProbability(G4int diQuark) const
{
  G4double sum = 0;
  typedef std::vector<G4SPPartonInfo *>::const_iterator iter;
  iter i;
  for (i = thePartonInfo.begin(); i!=thePartonInfo.end(); i++)
  {
    if (std::abs((*i)->GetDiQuark()) == std::abs(diQuark))
    {
      sum += (*i)->GetProbability();
    }
  }
  return sum;
}


G4int G4SPBaryon::
MatchDiQuarkAndGetQuark(const G4SPBaryon & aBaryon, G4int & aDiQuark) const
{
  G4int    result=0;
  typedef std::vector<G4SPPartonInfo *>::const_iterator iter;
  iter i;
  G4double running = 0;
  G4double total = 0;
  for (i = thePartonInfo.begin(); i!=thePartonInfo.end(); i++)
  {
    total += aBaryon.GetProbability((*i)->GetDiQuark());
  }
  G4double random = G4UniformRand();
  for(i = thePartonInfo.begin(); i!=thePartonInfo.end(); i++)
  {
    running += aBaryon.GetProbability((*i)->GetDiQuark());
    if (random<running/total)
    {
      result = (*i)->GetQuark();  // (diquark annihilated)
      aDiQuark = (*i)->GetDiQuark();
      break;
    }
  }
  return result;
}


void G4SPBaryon::
SampleQuarkAndDiquark(G4int & quark, G4int & diQuark) const
{
  typedef std::vector<G4SPPartonInfo *>::const_iterator iter;

  G4double random = G4UniformRand();
  G4double sum = 0;
  iter i;

  for (i=thePartonInfo.begin() ; i!=thePartonInfo.end(); i++)
  {
    sum += (*i)->GetProbability();
    if (sum > random)
    {
      if (theDefinition->GetPDGEncoding() < 0)
      {
        quark = (*i)->GetDiQuark();
        diQuark = (*i)->GetQuark();
      }
      else
      {
        quark = (*i)->GetQuark();
        diQuark = (*i)->GetDiQuark();
      }
      break;
    }
  }
}


void G4SPBaryon::
FindDiquark(G4int quark, G4int & diQuark) const
{
  typedef std::vector<G4SPPartonInfo *>::const_iterator iter;
  G4double sum = 0;
  iter i;
  for (i=thePartonInfo.begin() ; i!=thePartonInfo.end(); i++)
  {
    if (std::abs((*i)->GetQuark()) == std::abs(quark))
    {
      sum += (*i)->GetProbability();
    }
  }
  G4double random = G4UniformRand();
  G4double running = 0;
  for (i=thePartonInfo.begin() ; i!=thePartonInfo.end(); i++) {
    if (std::abs((*i)->GetQuark()) == std::abs(quark))
    {
      running += (*i)->GetProbability();
      if (running/sum >= random)
      {
        diQuark = (*i)->GetDiQuark();
        break;
      }
    }
  }
}


G4SPBaryon::
G4SPBaryon(G4Proton * aProton)
{
  theDefinition = aProton;
  thePartonInfo.push_back(new G4SPPartonInfo(2203, 1, 1./3./2.)); // uu_1, d 
  thePartonInfo.push_back(new G4SPPartonInfo(2103, 2, 1./6.*2.)); // ud_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(2101, 2, 1./2.));    // ud_0, u
}

G4SPBaryon::
G4SPBaryon(G4AntiProton * aAntiProton)
{
  theDefinition = aAntiProton;
  thePartonInfo.push_back(new G4SPPartonInfo(-2203, -1, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-2103, -2, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-2101, -2, 1./2.));
}


G4SPBaryon::
G4SPBaryon(G4Neutron * aNeutron)
{
  theDefinition = aNeutron;
  thePartonInfo.push_back(new G4SPPartonInfo(2103, 1, 1./6.*2.)); // ud_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(2101, 1, 1./2.   )); // ud_0, d
  thePartonInfo.push_back(new G4SPPartonInfo(1103, 2, 1./3./2 )); // dd_1, u
}

G4SPBaryon::
G4SPBaryon(G4AntiNeutron * aAntiNeutron)
{
  theDefinition = aAntiNeutron;
  thePartonInfo.push_back(new G4SPPartonInfo(-2103, -1, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-2101, -1, 1./2.));
  thePartonInfo.push_back(new G4SPPartonInfo(-1103, -2, 1./3.));
}


G4SPBaryon::
G4SPBaryon(G4Lambda * aLambda)
{
  theDefinition = aLambda;
  thePartonInfo.push_back(new G4SPPartonInfo(2103, 3, 1./3.));  // ud_1, s
  thePartonInfo.push_back(new G4SPPartonInfo(3203, 1, 1./4.));  // su_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(3201, 1, 1./12.)); // su_0, d
  thePartonInfo.push_back(new G4SPPartonInfo(3103, 2, 1./4.));  // sd_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(3101, 2, 1./12.)); // sd_0, u
}

G4SPBaryon::
G4SPBaryon(G4AntiLambda * aAntiLambda)
{
  theDefinition = aAntiLambda;
  thePartonInfo.push_back(new G4SPPartonInfo(-2103, -3, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3203, -1, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3201, -1, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3103, -2, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3101, -2, 1./12.));
}


G4SPBaryon::
G4SPBaryon(G4SigmaPlus * aSigmaPlus)
{
  theDefinition = aSigmaPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(2203, 3, 1./3.)); // uu_1, s
  thePartonInfo.push_back(new G4SPPartonInfo(3203, 2, 1./6.)); // su_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(3201, 2, 1./2.)); // su_0, u
}

G4SPBaryon::
G4SPBaryon(G4AntiSigmaPlus * aAntiSigmaPlus)
{
  theDefinition = aAntiSigmaPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(-2203, -3, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3203, -2, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3201, -2, 1./2.));
}


G4SPBaryon::
G4SPBaryon(G4SigmaZero * aSigmaZero)
{
  theDefinition = aSigmaZero;
  thePartonInfo.push_back(new G4SPPartonInfo(2103, 3, 1./3.));  // ud_1, s
  thePartonInfo.push_back(new G4SPPartonInfo(3203, 1, 1./12.)); // su_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(3201, 1, 1./4.));  // su_0, d
  thePartonInfo.push_back(new G4SPPartonInfo(3103, 2, 1./12.)); // sd_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(3101, 2, 1./4.));  // sd_0, u
}

G4SPBaryon::
G4SPBaryon(G4AntiSigmaZero * aAntiSigmaZero)
{
  theDefinition = aAntiSigmaZero;
  thePartonInfo.push_back(new G4SPPartonInfo(-2103, -3, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3203, -1, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3201, -1, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3103, -2, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3101, -2, 1./4.));
}


G4SPBaryon::
G4SPBaryon(G4SigmaMinus * aSigmaMinus)
{
  theDefinition = aSigmaMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(1103, 3, 1./3.)); // dd_1, s
  thePartonInfo.push_back(new G4SPPartonInfo(3103, 1, 1./6.)); // sd_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(3101, 1, 1./2.)); // sd_0, d
}

G4SPBaryon::
G4SPBaryon(G4AntiSigmaMinus * aAntiSigmaMinus)
{
  theDefinition = aAntiSigmaMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(-1103, -3, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3103, -1, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3101, -1, 1./2.));
}


G4SPBaryon::
G4SPBaryon(G4XiZero * aXiZero)
{
  theDefinition = aXiZero;
  thePartonInfo.push_back(new G4SPPartonInfo(3203, 3, 1./6.)); // su_1, s
  thePartonInfo.push_back(new G4SPPartonInfo(3201, 3, 1./2.)); // su_0, s
  thePartonInfo.push_back(new G4SPPartonInfo(3303, 2, 1./3.)); // ss_1, u
}

G4SPBaryon::
G4SPBaryon(G4AntiXiZero * aAntiXiZero)
{
  theDefinition = aAntiXiZero;
  thePartonInfo.push_back(new G4SPPartonInfo(-3203, -3, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3201, -3, 1./2.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3303, -2, 1./3.));
}


G4SPBaryon::
G4SPBaryon(G4XiMinus * aXiMinus)
{
  theDefinition = aXiMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(3103, 3, 1./6.)); // sd_1, s
  thePartonInfo.push_back(new G4SPPartonInfo(3101, 3, 1./2.)); // sd_0, s
  thePartonInfo.push_back(new G4SPPartonInfo(3303, 1, 1./3.)); // ss_1, d
}

G4SPBaryon::
G4SPBaryon(G4AntiXiMinus * aAntiXiMinus)
{
  theDefinition = aAntiXiMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(-3103, -3, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3101, -3, 1./2.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3303, -1, 1./3.));
}


G4SPBaryon::
G4SPBaryon(G4OmegaMinus * anOmegaMinus)
{
  theDefinition = anOmegaMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(3303, 3, 1.)); // ss_1, s
}

G4SPBaryon::
G4SPBaryon(G4AntiOmegaMinus * anAntiOmegaMinus)
{
  theDefinition = anAntiOmegaMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(-3303, -3, 1.));
}


// non static particles
G4SPBaryon::
G4SPBaryon(G4ParticleDefinition * aDefinition)
{
  theDefinition = aDefinition;
  if (theDefinition == G4ParticleTable::GetParticleTable()->FindParticle(2224)) // Delta++
  {
    thePartonInfo.push_back(new G4SPPartonInfo(2203, 2, 1.)); // uu_1, u
  }
  else if (theDefinition == G4ParticleTable::GetParticleTable()->FindParticle(-2224)) // anti Delta++
  {
    thePartonInfo.push_back(new G4SPPartonInfo(-2203, -2, 1.));
  }
  else if (theDefinition == G4ParticleTable::GetParticleTable()->FindParticle(2214)) // Delta+
  {
    thePartonInfo.push_back(new G4SPPartonInfo(2203, 1, 1./3.)); // uu_1, d
    thePartonInfo.push_back(new G4SPPartonInfo(2103, 2, 2./3.)); // ud_1, u
  }
  else if (theDefinition == G4ParticleTable::GetParticleTable()->FindParticle(-2214)) // anti Delta+
  {
    thePartonInfo.push_back(new G4SPPartonInfo(-2203, -1, 1./3.));
    thePartonInfo.push_back(new G4SPPartonInfo(-2103, -2, 2./3.));
  }
  else if (theDefinition == G4ParticleTable::GetParticleTable()->FindParticle(2114)) // Delta0
  {
    thePartonInfo.push_back(new G4SPPartonInfo(2103, 1, 2./3.)); // ud_1, d
    thePartonInfo.push_back(new G4SPPartonInfo(1103, 2, 1./3.)); // dd_1, u
  }
  else if (theDefinition == G4ParticleTable::GetParticleTable()->FindParticle(-2114)) // anti Delta0
  {
    thePartonInfo.push_back(new G4SPPartonInfo(-2103, -1, 2./3.));
    thePartonInfo.push_back(new G4SPPartonInfo(-2103, -2, 1./3.));
  }
  else if (theDefinition == G4ParticleTable::GetParticleTable()->FindParticle(1114)) // Delta-
  {
    thePartonInfo.push_back(new G4SPPartonInfo(1103, 1, 1.)); // dd_1, d
  }
  else if (theDefinition == G4ParticleTable::GetParticleTable()->FindParticle(-1114)) // anti Delta-
  {
    thePartonInfo.push_back(new G4SPPartonInfo(-1103, -1, 1.));
  }  
}


G4SPBaryon::~G4SPBaryon()
{
  for (unsigned int i=0;i<thePartonInfo.size(); i++) delete thePartonInfo[i];
}


// Extension to charmed and bottom baryons and anti-baryons
//    G4SPPartonInfo(G4int diq, G4int q, G4double prob)

G4SPBaryon::G4SPBaryon(G4LambdacPlus * aLambdacPlus) {
  // lambda_c+(udc) treated as lambda(uds) with s replaced by c.
  theDefinition = aLambdacPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(2103, 4, 1./3.));  // ud_1, c
  thePartonInfo.push_back(new G4SPPartonInfo(4203, 1, 1./4.));  // cu_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(4201, 1, 1./12.)); // cu_0, d
  thePartonInfo.push_back(new G4SPPartonInfo(4103, 2, 1./4.));  // cd_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(4101, 2, 1./12.)); // cd_0, u
}

G4SPBaryon::G4SPBaryon(G4AntiLambdacPlus * aAntiLambdacPlus) {
  theDefinition = aAntiLambdacPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(-2103, -4, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4203, -1, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4201, -1, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4103, -2, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4101, -2, 1./12.));
}


G4SPBaryon::G4SPBaryon(G4SigmacPlusPlus * aSigmacPlusPlus) {
  // sigma_c++(uuc) treated as sigma+(uus) with s replaced by c.
  theDefinition = aSigmacPlusPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(2203, 4, 1./3.)); // uu_1, c
  thePartonInfo.push_back(new G4SPPartonInfo(4203, 2, 1./6.)); // cu_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(4201, 2, 1./2.)); // cu_0, u
}

G4SPBaryon::G4SPBaryon(G4AntiSigmacPlusPlus * aAntiSigmacPlusPlus) {
  theDefinition = aAntiSigmacPlusPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(-2203, -4, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4203, -2, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4201, -2, 1./2.));
}


G4SPBaryon::G4SPBaryon(G4SigmacPlus * aSigmacPlus) {
  // sigma_c+(udc) treated as sigma0(uds) with s replaced by c.
  theDefinition = aSigmacPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(2103, 4, 1./3.));  // ud_1, c
  thePartonInfo.push_back(new G4SPPartonInfo(4203, 1, 1./12.)); // cu_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(4201, 1, 1./4.));  // cu_0, d
  thePartonInfo.push_back(new G4SPPartonInfo(4103, 2, 1./12.)); // cd_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(4101, 2, 1./4.));  // cd_0, u
}

G4SPBaryon::G4SPBaryon(G4AntiSigmacPlus * aAntiSigmacPlus) {
  theDefinition = aAntiSigmacPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(-2103, -4, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4203, -1, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4201, -1, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4103, -2, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4101, -2, 1./4.));
}


G4SPBaryon::G4SPBaryon(G4SigmacZero * aSigmacZero) {
  // sigma_c0(ddc) treated as sigma-(dds) replacing s with c.
  theDefinition = aSigmacZero;
  thePartonInfo.push_back(new G4SPPartonInfo(1103, 4, 1./3.)); // dd_1, c
  thePartonInfo.push_back(new G4SPPartonInfo(4103, 1, 1./6.)); // cd_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(4101, 1, 1./2.)); // cd_0, d
}

G4SPBaryon::G4SPBaryon(G4AntiSigmacZero * aAntiSigmacZero) {
  theDefinition = aAntiSigmacZero;
  thePartonInfo.push_back(new G4SPPartonInfo(-1103, -4, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4103, -1, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4101, -1, 1./2.));
}


G4SPBaryon::G4SPBaryon(G4XicPlus * aXicPlus) {
  // xi_c+(usc) treated as xi0(uss) replacing s with c.
  theDefinition = aXicPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(3203, 4, 1./6.)); // su_1, c
  thePartonInfo.push_back(new G4SPPartonInfo(3201, 4, 1./2.)); // su_0, c
  thePartonInfo.push_back(new G4SPPartonInfo(4303, 2, 1./3.)); // cs_1, u  
}

G4SPBaryon::G4SPBaryon(G4AntiXicPlus * aAntiXicPlus) {
  theDefinition = aAntiXicPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(-3203, -4, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3201, -4, 1./2.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4303, -2, 1./3.));
}


G4SPBaryon::G4SPBaryon(G4XicZero * aXicZero) {
  // xi_c0(dsc) treated as xi-(dss) replacing s with c.
  theDefinition = aXicZero;
  thePartonInfo.push_back(new G4SPPartonInfo(3103, 4, 1./6.)); // sd_1, c
  thePartonInfo.push_back(new G4SPPartonInfo(3101, 4, 1./2.)); // sd_0, c
  thePartonInfo.push_back(new G4SPPartonInfo(4303, 1, 1./3.)); // cs_1, d
}

G4SPBaryon::G4SPBaryon(G4AntiXicZero * aAntiXicZero) {
  theDefinition = aAntiXicZero;
  thePartonInfo.push_back(new G4SPPartonInfo(-3103, -4, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3101, -4, 1./2.));
  thePartonInfo.push_back(new G4SPPartonInfo(-4303, -1, 1./3.));
}


G4SPBaryon::G4SPBaryon(G4OmegacZero * aOmegacZero) {
  // omega_c0(ssc) treated as omega-(sss) with s replaced by c.
  theDefinition = aOmegacZero;
  thePartonInfo.push_back(new G4SPPartonInfo(3303, 4, 1.)); // ss_1, c
}

G4SPBaryon::G4SPBaryon(G4AntiOmegacZero * aAntiOmegacZero) {
  theDefinition = aAntiOmegacZero;
  thePartonInfo.push_back(new G4SPPartonInfo(-3303, -4, 1.));
}


G4SPBaryon::G4SPBaryon(G4Lambdab * aLambdab) {
  // lambda_b(udb) treated as lambda-(uds) replacing s with b.
  theDefinition = aLambdab;
  thePartonInfo.push_back(new G4SPPartonInfo(2103, 5, 1./3.));  // ud_1, b
  thePartonInfo.push_back(new G4SPPartonInfo(5203, 1, 1./4.));  // bu_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(5201, 1, 1./12.)); // bu_0, d
  thePartonInfo.push_back(new G4SPPartonInfo(5103, 2, 1./4.));  // bd_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(5101, 2, 1./12.)); // bd_0, u
}

G4SPBaryon::G4SPBaryon(G4AntiLambdab * aAntiLambdab) {
  theDefinition = aAntiLambdab;
  thePartonInfo.push_back(new G4SPPartonInfo(-2103, -5, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5203, -1, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5201, -1, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5103, -2, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5101, -2, 1./12.));
}


G4SPBaryon::G4SPBaryon(G4SigmabPlus * aSigmabPlus) {
  // sigma_b+(uub) treated as sigma+(uus) replacing s with b.
  theDefinition = aSigmabPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(2203, 5, 1./3.)); // uu_1, b
  thePartonInfo.push_back(new G4SPPartonInfo(5203, 2, 1./6.)); // bu_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(5201, 2, 1./2.)); // bu_0, u
}

G4SPBaryon::G4SPBaryon(G4AntiSigmabPlus * aAntiSigmabPlus) {
  theDefinition = aAntiSigmabPlus;
  thePartonInfo.push_back(new G4SPPartonInfo(-2203, -5, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5203, -2, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5201, -2, 1./2.));
}


G4SPBaryon::G4SPBaryon(G4SigmabZero * aSigmabZero) {
  // sigma_b0(udb) treated as sigma0(uds) replacing s with b.
  theDefinition = aSigmabZero;
  thePartonInfo.push_back(new G4SPPartonInfo(2103, 5, 1./3.));  // ud_1, b
  thePartonInfo.push_back(new G4SPPartonInfo(5203, 1, 1./12.)); // bu_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(5201, 1, 1./4.));  // bu_0, d
  thePartonInfo.push_back(new G4SPPartonInfo(5103, 2, 1./12.)); // bd_1, u
  thePartonInfo.push_back(new G4SPPartonInfo(5101, 2, 1./4.));  // bd_0, u
}

G4SPBaryon::G4SPBaryon(G4AntiSigmabZero * aAntiSigmabZero) {
  theDefinition = aAntiSigmabZero;
  thePartonInfo.push_back(new G4SPPartonInfo(-2103, -5, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5203, -1, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5201, -1, 1./4.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5103, -2, 1./12.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5101, -2, 1./4.));
}


G4SPBaryon::G4SPBaryon(G4SigmabMinus * aSigmabMinus) {
  // sigma_b-(ddb) treated as sigma-(dds) replacing s with b.
  theDefinition = aSigmabMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(1103, 5, 1./3.)); // dd_1, b
  thePartonInfo.push_back(new G4SPPartonInfo(5103, 1, 1./6.)); // bd_1, d
  thePartonInfo.push_back(new G4SPPartonInfo(5101, 1, 1./2.)); // bd_0, d
}

G4SPBaryon::G4SPBaryon(G4AntiSigmabMinus * aAntiSigmabMinus) {
  theDefinition = aAntiSigmabMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(-1103, -5, 1./3.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5103, -1, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5101, -1, 1./2.));
}


G4SPBaryon::G4SPBaryon(G4XibZero * aXibZero) {
  // xi_b0(usb) treated as xi0(uss) replacing s with b.
  theDefinition = aXibZero;
  thePartonInfo.push_back(new G4SPPartonInfo(3203, 5, 1./6.)); // su_1, b
  thePartonInfo.push_back(new G4SPPartonInfo(3201, 5, 1./2.)); // su_0, b
  thePartonInfo.push_back(new G4SPPartonInfo(5303, 2, 1./3.)); // bs_1, u
}

G4SPBaryon::G4SPBaryon(G4AntiXibZero * aAntiXibZero) {
  theDefinition = aAntiXibZero;
  thePartonInfo.push_back(new G4SPPartonInfo(-3203, -5, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3201, -5, 1./2.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5303, -2, 1./3.));
}


G4SPBaryon::G4SPBaryon(G4XibMinus * aXibMinus) {
  // xi_b-(dsb) treated as xi-(dss) replacing s with b.
  theDefinition = aXibMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(3103, 5, 1./6.)); // sd_1, b
  thePartonInfo.push_back(new G4SPPartonInfo(3101, 5, 1./2.)); // sd_0, b
  thePartonInfo.push_back(new G4SPPartonInfo(5303, 1, 1./3.)); // bs_1, d
}

G4SPBaryon::G4SPBaryon(G4AntiXibMinus * aAntiXibMinus) {
  theDefinition = aAntiXibMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(-3103, -5, 1./6.));
  thePartonInfo.push_back(new G4SPPartonInfo(-3101, -5, 1./2.));
  thePartonInfo.push_back(new G4SPPartonInfo(-5303, -1, 1./3.));
}


G4SPBaryon::G4SPBaryon(G4OmegabMinus * aOmegabMinus) {
  // omega_b-(ssb) treated as omega-(sss) replacing s with b.
  theDefinition = aOmegabMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(3303, 5, 1.)); // ss_1, b
}

G4SPBaryon::G4SPBaryon(G4AntiOmegabMinus * aAntiOmegabMinus) {
  theDefinition = aAntiOmegabMinus;
  thePartonInfo.push_back(new G4SPPartonInfo(-3303, -5, 1.));
}
