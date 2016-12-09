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
#ifndef G4PomeronCrossSection_h
#define G4PomeronCrossSection_h 1
//
// $Id: G4PomeronCrossSection.hh 100828 2016-11-02 15:25:59Z gcosmo $
//
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4Gamma.hh"

class G4PomeronCrossSection
{
  public:
    G4PomeronCrossSection(const G4ParticleDefinition * );
    G4PomeronCrossSection(const G4Proton * );
    G4PomeronCrossSection(const G4Neutron * );

    G4PomeronCrossSection(const G4PionPlus * );
    G4PomeronCrossSection(const G4PionMinus * );
    G4PomeronCrossSection(const G4PionZero * );

    G4PomeronCrossSection(const G4KaonPlus * );
    G4PomeronCrossSection(const G4KaonMinus * );
    G4PomeronCrossSection(const G4KaonZero * );
    G4PomeronCrossSection(const G4KaonZeroLong * );
    G4PomeronCrossSection(const G4KaonZeroShort * );

    G4PomeronCrossSection(const G4Gamma * );

    ~G4PomeronCrossSection();

    // s = (center of mass energy)**2	
    G4double GetTotalCrossSection(const G4double s);
    G4double GetDiffractiveCrossSection(const G4double s);
    G4double GetElasticCrossSection(const G4double s);
    G4double GetInelasticCrossSection(const G4double s);

    G4double GetTotalProbability(const G4double s, const G4double impactsquare);
    G4double GetDiffractiveProbability(const G4double s, const G4double impactsquare);
    G4double GetNondiffractiveProbability(const G4double s, const G4double impactsquare);
    G4double GetElasticProbability(const G4double s, const G4double impactsquare);
    G4double GetInelasticProbability(const G4double s, const G4double impactsquare);
				     
    G4double GetCutPomeronProbability(const G4double s, const G4double impactsquare, 
                                      const G4int nPomerons);
	
    void Setgamma(const G4double agam); // temporary only! GF.
    G4double SoftEikonal(G4double s, G4double impactsquare);
    G4double HardEikonal(G4double s, G4double impactsquare);
        
    void Pomeron_S(G4double apomeron_S){ pomeron_S = apomeron_S;}
    void Pomeron_Gamma(G4double apomeron_Gamma){ pomeron_Gamma = apomeron_Gamma;}
    void Pomeron_C(G4double apomeron_C){ pomeron_C = apomeron_C;}
    void Pomeron_Rsquare(G4double apomeron_Rsquare){ pomeron_Rsquare = apomeron_Rsquare;}
    void Pomeron_Alpha(G4double apomeron_Alpha){ pomeron_Alpha = apomeron_Alpha;}
    void Pomeron_Alphaprime(G4double apomeron_Alphaprime){ pomeron_Alphaprime = apomeron_Alphaprime;}
    void Pomeron_Gamma_Hard(G4double apomeron_Gamma_Hard){ pomeron_Gamma_Hard = apomeron_Gamma_Hard;}
    void Pomeron_Alpha_Hard(G4double apomeron_Alpha_Hard){ pomeron_Alpha_Hard = apomeron_Alpha_Hard;}

  private: 
    G4double PowerSoft(const G4double s);
    G4double PowerHard(const G4double s);
    G4double LambdaSoft(const G4double s);
    G4double LambdaHard(const G4double s);
    G4double Zsoft(const G4double s);
    G4double Zhard(const G4double s);
  
    G4PomeronCrossSection();
    void InitForNucleon();
    void InitForPion();
    void InitForKaon();
    void InitForGamma();
	
    G4double Expand(G4double z);
    inline G4double Z(const G4double Scms);
    inline G4double SigP(const G4double Scms);
    inline G4double Power(const G4double Scms);
    inline G4double Lambda(const G4double s);
    inline G4double Eikonal(const G4double s,const G4double impactsquare);
        
    G4double pomeron_Alpha;
    G4double pomeron_Alpha_Hard;
    G4double pomeron_Alphaprime;
    G4double pomeron_C;
    G4double pomeron_Gamma;
    G4double pomeron_Gamma_Hard;
    G4double pomeron_Rsquare;
    G4double pomeron_S;
};

#endif

