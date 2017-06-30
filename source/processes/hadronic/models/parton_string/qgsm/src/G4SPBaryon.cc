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
	for(i = thePartonInfo.begin(); i!=thePartonInfo.end(); i++)
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
	for(i = thePartonInfo.begin(); i!=thePartonInfo.end(); i++)
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
	for(i = thePartonInfo.begin(); i!=thePartonInfo.end(); i++)
	{
		total += aBaryon.GetProbability((*i)->GetDiQuark());
	}
	G4double random = G4UniformRand();                            //*total; Vova 17 Sept.
	for(i = thePartonInfo.begin(); i!=thePartonInfo.end(); i++)
	{
		running += aBaryon.GetProbability((*i)->GetDiQuark());
		if(random<running/total)                             // if(random/total<running) Vova 17 Sept.
		{
			result = (*i)->GetQuark(); // (diquark annihilated)
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

	for(i=thePartonInfo.begin() ; i!=thePartonInfo.end(); i++)
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
	for(i=thePartonInfo.begin() ; i!=thePartonInfo.end(); i++)
	{
		if (std::abs((*i)->GetQuark()) == std::abs(quark))
		{
			sum += (*i)->GetProbability();
		}
	}
	G4double random = G4UniformRand();
	G4double running = 0;
	for(i=thePartonInfo.begin() ; i!=thePartonInfo.end(); i++)
	{
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
	theDefinition = aProton;                              // Uzhi
//
	thePartonInfo.push_back(new G4SPPartonInfo(2203, 1, 1./3./2.)); // uu_1, d 
	thePartonInfo.push_back(new G4SPPartonInfo(2103, 2, 1./6.*2.)); // ud_1, u
	thePartonInfo.push_back(new G4SPPartonInfo(2101, 2, 1./2.));    // ud_0, u
//
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
	theDefinition = aNeutron;                                // Uzhi
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
	thePartonInfo.push_back(new G4SPPartonInfo(2203, 3, 1./3.));
	thePartonInfo.push_back(new G4SPPartonInfo(3203, 2, 1./6.));
	thePartonInfo.push_back(new G4SPPartonInfo(3201, 2, 1./2.));
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
	thePartonInfo.push_back(new G4SPPartonInfo(2103, 3, 1./3.));
	thePartonInfo.push_back(new G4SPPartonInfo(3203, 1, 1./12.));
	thePartonInfo.push_back(new G4SPPartonInfo(3201, 1, 1./4.));
	thePartonInfo.push_back(new G4SPPartonInfo(3103, 2, 1./12.));
	thePartonInfo.push_back(new G4SPPartonInfo(3101, 2, 1./4.));
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
	thePartonInfo.push_back(new G4SPPartonInfo(1103, 3, 1./3.));
	thePartonInfo.push_back(new G4SPPartonInfo(3103, 1, 1./6.));
	thePartonInfo.push_back(new G4SPPartonInfo(3101, 1, 1./2.));
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
G4SPBaryon(G4XiMinus * aXiMinus)
{
	theDefinition = aXiMinus;
	thePartonInfo.push_back(new G4SPPartonInfo(3103, 3, 1./6.));
	thePartonInfo.push_back(new G4SPPartonInfo(3101, 3, 1./2.));
	thePartonInfo.push_back(new G4SPPartonInfo(3303, 1, 1./3.));
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
G4SPBaryon(G4XiZero * aXiZero)
{
	theDefinition = aXiZero;
	thePartonInfo.push_back(new G4SPPartonInfo(3203, 3, 1./6.));
	thePartonInfo.push_back(new G4SPPartonInfo(3201, 3, 1./2.));
	thePartonInfo.push_back(new G4SPPartonInfo(3303, 2, 1./3.));
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
G4SPBaryon(G4OmegaMinus * anOmegaMinus)
{
	theDefinition = anOmegaMinus;
	thePartonInfo.push_back(new G4SPPartonInfo(3303, 3, 1.));
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
	if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(2224))// D++
	{
		thePartonInfo.push_back(new G4SPPartonInfo(2203, 2, 1.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-2224))// anti D++
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-2203, -2, 1.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(2214))// D+
	{
		thePartonInfo.push_back(new G4SPPartonInfo(2203, 1, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(2103, 2, 2./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-2214))// anti D+
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-2203, -1, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(-2103, -2, 2./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(2114))// D0
	{
		thePartonInfo.push_back(new G4SPPartonInfo(2103, 1, 2./3.));
// Uzhi		thePartonInfo.push_back(new G4SPPartonInfo(2103, 2, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(1103, 2, 1./3.));  // Uzhi 14.05.2014
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-2114))// anti D0
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-2103, -1, 2./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(-2103, -2, 1./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(1114))// D-
	{
		thePartonInfo.push_back(new G4SPPartonInfo(1103, 1, 1.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-1114))// anti D-
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-1103, -1, 1.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(3224))// S*+
	{
		thePartonInfo.push_back(new G4SPPartonInfo(2203, 3, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(3203, 2, 2./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-3224))// anti S*+
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-2203, -3, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(-3203, -2, 2./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(3214))// S*0
	{
		thePartonInfo.push_back(new G4SPPartonInfo(2103, 3, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(3203, 1, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(3103, 2, 1./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-3214))// anti S*0
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-2103, -3, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(-3203, -1, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(-3103, -2, 1./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(3114))// S*-
	{
		thePartonInfo.push_back(new G4SPPartonInfo(1103, 3, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(3103, 1, 2./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-3224))// anti S*-
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-1103, -3, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(-3103, -1, 2./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(3324))// Xi*0
	{
		thePartonInfo.push_back(new G4SPPartonInfo(3203, 3, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(3303, 2, 2./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-3324))// anti Xi*0
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-3203, -3, 1./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(-3303, -2, 2./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(3314))// Xi*-
	{
		thePartonInfo.push_back(new G4SPPartonInfo(3103, 3, 2./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(3303, 1, 1./3.));
	}
	else if(theDefinition ==
			G4ParticleTable::GetParticleTable()->FindParticle(-3314))// anti Xi*-
	{
		thePartonInfo.push_back(new G4SPPartonInfo(-3103, -3, 2./3.));
		thePartonInfo.push_back(new G4SPPartonInfo(-3303, -1, 1./3.));
	}
}

G4SPBaryon::~G4SPBaryon()
{
	for(unsigned int i=0;i<thePartonInfo.size(); i++) delete thePartonInfo[i];
}
