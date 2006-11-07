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
//
// $Id: G4PomeronCrossSection.cc,v 1.6 2006-11-07 12:51:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4PomeronCrossSection.hh"

G4PomeronCrossSection::G4PomeronCrossSection()
{;}


G4PomeronCrossSection::~G4PomeronCrossSection()
{;}
//**********************************************************************************************

G4PomeronCrossSection::G4PomeronCrossSection(const G4ParticleDefinition * particle)
    {
    G4int Encoding = std::abs(particle->GetPDGEncoding());

    if (std::abs(particle->GetBaryonNumber())!=0)
	InitForNucleon();
    else if (Encoding/100== 3 || Encoding/10 == 3)
	InitForKaon();
    else
	InitForPion();
    }

//**********************************************************************************************

G4PomeronCrossSection::G4PomeronCrossSection(const G4Proton * )
{
	InitForNucleon();
}
G4PomeronCrossSection::G4PomeronCrossSection(const G4Neutron * )
{
	InitForNucleon();
}

G4PomeronCrossSection::G4PomeronCrossSection(const G4PionPlus * )
{
	InitForPion();
}
G4PomeronCrossSection::G4PomeronCrossSection(const G4PionMinus * )
{
	InitForPion();
}
G4PomeronCrossSection::G4PomeronCrossSection(const G4PionZero * )
{
	InitForPion();
}

G4PomeronCrossSection::G4PomeronCrossSection(const G4KaonPlus * )
{
	InitForKaon();
}
G4PomeronCrossSection::G4PomeronCrossSection(const G4KaonMinus * )
{
	InitForKaon();
}
G4PomeronCrossSection::G4PomeronCrossSection(const G4KaonZero * )
{
	InitForKaon();
}
G4PomeronCrossSection::G4PomeronCrossSection(const G4KaonZeroLong * )
{
	InitForKaon();
}
G4PomeronCrossSection::G4PomeronCrossSection(const G4KaonZeroShort * )
{
	InitForKaon();
}
G4PomeronCrossSection::G4PomeronCrossSection(const G4Gamma * )
{
	InitForGamma();
}

G4double G4PomeronCrossSection::GetTotalCrossSection(const G4double s)
{
	G4double FZ2= Expand(Z(s)/2);
	return SigP(s) * FZ2;
}


G4double G4PomeronCrossSection::GetElasticCrossSection(const G4double s)
{
	return SigP(s)/pomeron_C *(Expand(Z(s)/2) - Expand(Z(s)));
}

G4double G4PomeronCrossSection::GetDiffractiveCrossSection(const G4double s)
{
	return ( pomeron_C -1) * GetElasticCrossSection(s);
}

G4double G4PomeronCrossSection::GetInelasticCrossSection(const G4double s)
{
	return GetTotalCrossSection(s) - GetElasticCrossSection(s);
}

//-------------------------Probabilities ----------------------------

G4double G4PomeronCrossSection::GetTotalProbability(const G4double s, 
				  const G4double impactsquare)
{
	return 2/pomeron_C*(1-std::exp(-1*Eikonal(s,impactsquare)));	
}

G4double G4PomeronCrossSection::GetDiffractiveProbability(const G4double s, 
				  const G4double impactsquare)
{
	return (pomeron_C-1)/pomeron_C * 
	                     (GetTotalProbability(s,impactsquare) -
			      GetNondiffractiveProbability(s,impactsquare));
}

G4double G4PomeronCrossSection::GetNondiffractiveProbability(const G4double s, 
				  const G4double impactsquare)
{
	return (1-std::exp(-2*Eikonal(s,impactsquare)))/pomeron_C;	
}

G4double G4PomeronCrossSection::GetElasticProbability(const G4double s, 
				  const G4double impactsquare)
{
	return (GetTotalProbability(s,impactsquare) -
	        GetInelasticProbability(s,impactsquare));
}

G4double G4PomeronCrossSection::GetInelasticProbability(const G4double s, 
				  const G4double impactsquare)
{
	return GetNondiffractiveProbability(s,impactsquare) +
		GetDiffractiveProbability(s,impactsquare);
}

G4double G4PomeronCrossSection::GetCutPomeronProbability(const G4double s,
			    const G4double impactsquare, const G4int nPomerons)
{
	G4double factorial=1;
	for (G4int i=2; i<= nPomerons ; i++) 
	{   factorial *= i;
	}
	 
	return std::exp(-2*Eikonal(s,impactsquare))/pomeron_C*
	        std::pow(2*Eikonal(s,impactsquare),nPomerons)/factorial;
}

// ---------------Temporary --- GF
void G4PomeronCrossSection::Setgamma(const G4double agam)
{
	pomeron_Gamma=agam/GeV/GeV;
}


//-----------------  private/Implementation methods

void G4PomeronCrossSection::InitForNucleon()
{
//	pomeron_S=		3.0*GeV*GeV;
	pomeron_S=		2.7*GeV*GeV;
//	pomeron_Gamma=		2.16/GeV/GeV;
//	pomeron_Gamma=		3.96/GeV/GeV;
	pomeron_Gamma=		(2.6+3.96)/GeV/GeV;
	pomeron_C=		1.4;
	pomeron_Rsquare=	3.56/GeV/GeV;
//	pomeron_Alpha=		1.0808;
	pomeron_Alpha=		0.9808;
	pomeron_Alphaprime=	0.25/GeV/GeV;     
        pomeron_Gamma_Hard =    0.0002/GeV/GeV;  // Note! if pomeron_Gamma_Hard != 0 to fit total pp-crosscection
                                                 // pomeron_Gamma_Soft shold be 2.35/GeV/GeV
        pomeron_Alpha_Hard =    1.47;            
}

void G4PomeronCrossSection::InitForPion()
{
	pomeron_S=		1.5*GeV*GeV;
//	pomeron_Gamma=		1.46/GeV/GeV;
	pomeron_Gamma=		2.17/GeV/GeV;
	pomeron_C=		1.6;
	pomeron_Rsquare=	2.36/GeV/GeV;
	pomeron_Alpha=		1.0808;
	pomeron_Alphaprime=	0.25/GeV/GeV;
        pomeron_Gamma_Hard =    0.0002/GeV/GeV;
        pomeron_Alpha_Hard =    1.47;
}

void G4PomeronCrossSection::InitForKaon()
{
	pomeron_S=		2.3*GeV*GeV;
//	pomeron_Gamma=		1.31/GeV/GeV;
	pomeron_Gamma=		1.92/GeV/GeV;
	pomeron_C=		1.8;
	pomeron_Rsquare=	1.96/GeV/GeV;
	pomeron_Alpha=		1.0808;
	pomeron_Alphaprime=	0.25/GeV/GeV;
        pomeron_Gamma_Hard =    0.0002/GeV/GeV;
        pomeron_Alpha_Hard =    1.47;
}

void G4PomeronCrossSection::InitForGamma()
{
	pomeron_S=		1.7*GeV*GeV;
//	pomeron_Gamma=		1.42/GeV/GeV;
	pomeron_Gamma=		2.07/GeV/GeV;
	pomeron_C=		1.7;
	pomeron_Rsquare=	2.16/GeV/GeV;
	pomeron_Alpha=		1.0808;
	pomeron_Alphaprime=	0.25/GeV/GeV;
        pomeron_Gamma_Hard =    0.0002/GeV/GeV;
        pomeron_Alpha_Hard =    1.47;
}

G4double G4PomeronCrossSection::Expand(G4double z)
{

	G4double sum=1.;
	G4double current=1.;
	for (G4int j=2; j<21; j++ )
	{
		current *= -z *(j-1)/sqr(j);
		sum+=current;
	}
	return sum;
}

inline G4double G4PomeronCrossSection::Power(const G4double s)
{
	return pomeron_Gamma *std::pow(s/pomeron_S, pomeron_Alpha -1);
}

inline G4double G4PomeronCrossSection::Z(const G4double s)
{
	return 2*pomeron_C * Power(s) / Lambda(s);
}

inline G4double G4PomeronCrossSection::Lambda(const G4double s)
{
	return pomeron_Rsquare+pomeron_Alphaprime*std::log(s/pomeron_S);
}

inline G4double G4PomeronCrossSection::SigP(const G4double s)
{
	return 8 * pi * hbarc_squared * Power(s);
}

inline G4double G4PomeronCrossSection::Eikonal(const G4double s,
					       const G4double impactsquare)
{
	return Z(s)/2 * std::exp(-impactsquare/(4*Lambda(s)*hbarc_squared));
}
//*************************************************************************************************
inline G4double G4PomeronCrossSection::PowerSoft(const G4double s)
    {
    return pomeron_Gamma *std::pow(s/pomeron_S, pomeron_Alpha -1);
    }

inline G4double G4PomeronCrossSection::PowerHard(const G4double s)
    {
    return pomeron_Gamma_Hard*std::pow(s/pomeron_S, pomeron_Alpha_Hard -1);
    }

inline G4double G4PomeronCrossSection::LambdaSoft(const G4double s)
    {
    return pomeron_Rsquare+pomeron_Alphaprime*std::log(s/pomeron_S);
    }
    
inline G4double G4PomeronCrossSection::LambdaHard(const G4double /*s*/)
    {
    return pomeron_Rsquare; //+pomeron_Alphaprime*std::log(s/pomeron_S);
    }

inline G4double G4PomeronCrossSection::Zsoft(const G4double s)
   {
   return 2*pomeron_C*PowerHard(s) / LambdaSoft(s);
   }

inline G4double G4PomeronCrossSection::Zhard(const G4double s)
   {
   return 2*pomeron_C*PowerHard(s)/LambdaHard(s);
   }

G4double G4PomeronCrossSection::SoftEikonal(G4double s, G4double impactsquare)
    {
    return Zsoft(s)/2*std::exp(-impactsquare/LambdaSoft(s)/hbarc_squared/4);
    }
     
G4double G4PomeronCrossSection::HardEikonal(G4double s, G4double impactsquare)
    {
    return Zhard(s)/2*std::exp(-impactsquare/LambdaHard(s)/hbarc_squared/4);
    }
    
//*************************************************************************************************
