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
//	G4ScreeningMottCrossSection.hh
//-------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4ScreeningMottCrossSection
//
// Author:      Cristina Consolandi
//
// Creation date: 20.10.2011
//
// Modifications:
// 27-05-2012 Added Analytic Fitting to the Mott Cross Section by means of G4MottCoefficients class.
//
//
// Class Description:
//      Computation of electron Coulomb Scattering Cross Section.
//      Suitable for high energy electrons and light target materials.
//
//      Reference:
//      M.J. Boschini et al.
//     "Non Ionizing Energy Loss induced by Electrons in the Space Environment"
//      Proc. of the 13th International Conference on Particle Physics and Advanced Technology
//      (13th ICPPAT, Como 3-7/10/2011), World Scientific (Singapore).
//      Available at: http://arxiv.org/abs/1111.4042v4
//
//      1) Mott Differential Cross Section Approximation:
//         For Target material up to Z=92 (U):
//         As described in http://arxiv.org/abs/1111.4042v4
//         par. 2.1 , eq. (16)-(17)
//         Else (Z>92):
//         W. A. McKinley and H. Fashbach, Phys. Rev. 74, (1948) 1759.
//      2) Screening coefficient:
//      vomn G. Moliere, Z. Naturforsh A2 (1947), 133-145; A3 (1948), 78.
//      3) Nuclear Form Factor:
//      A.V. Butkevich et al. Nucl. Instr. and Meth. in Phys. Res. A 488 (2002), 282-294.
//
// ----------------------------------------------------------------------------------------

//
#ifndef G4ScreeningMottCrossSection_h
#define G4ScreeningMottCrossSection_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "G4MottCoefficients.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4NistManager.hh"
#include "G4ThreeVector.hh"
#include "G4Pow.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static const G4int DIM = 750;

class G4ParticleDefinition;

class G4ScreeningMottCrossSection
{

public:

  	explicit G4ScreeningMottCrossSection();

  	virtual ~G4ScreeningMottCrossSection();

   	void Initialise(const G4ParticleDefinition*, G4double cosThetaLim);

	G4double GetScreeningAngle();
	void SetScreeningCoefficient();

  	inline void SetupParticle(const G4ParticleDefinition*);
 	void SetupKinematic(G4double kinEnergy ,G4double Z);

        G4double NuclearCrossSection(G4int, G4int);
        G4ThreeVector GetNewDirection();

	inline G4double GetMom2CM()const;
	inline G4double GetMom2Lab()const;
	inline G4double GetTrec() const;
        inline G4double GetScreeningCoefficient() const;
        inline G4double GetTotalCross() const;

        G4double McFcorrection(G4double);
	G4double RatioMottRutherford(G4double);
        G4double FormFactor2ExpHof(G4double);
        G4double FormFactor2Gauss(G4double);
	G4double FormFactor2UniformHelm(G4double);
	G4double GetScatteringAngle(G4int, G4int);
        G4double SetDifferentialXSection(G4double, G4double, G4int);

private:

 	G4ScreeningMottCrossSection & operator=(const  G4ScreeningMottCrossSection &right) = delete;
  	G4ScreeningMottCrossSection(const  G4ScreeningMottCrossSection&) = delete;

  	G4NistManager*  fNistManager;
	G4Pow*          fG4pow;
	G4MottCoefficients * mottcoeff;

	G4double 	        TotalCross;
	//cost - min - max
        G4double                cosThetaMin;// def 1.0
        G4double                cosThetaMax;// def -1.0

	G4double      	        cosTetMinNuc;
        G4double	        cosTetMaxNuc;

	//energy cut
  	G4double                ecut;
   	G4double                etag;

	// projectile........................
 	const G4ParticleDefinition* particle;

 	G4double                spin;
  	G4double                mass;

	//lab of incedent particle
        G4double                tkinLab;
        G4double                momLab2;
        G4double                invbetaLab2;

	//relative system with nucleus
	G4double 		mu_rel;
 	G4double                tkin;
  	G4double                mom2;
  	G4double                invbeta2;
	G4double		beta;
	G4double		gamma;

  	// target nucleus
  	G4double                targetZ;
	G4double 		targetA;
  	G4double                targetMass;
	G4double 		Trec;
 	G4double 		As;
	G4double 		coeffb[5][6];

	//constants
	G4double                 alpha;
  	G4double                 htc2;
        G4double                 e2;

        //angle
        G4double cross[DIM];
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4ScreeningMottCrossSection::SetupParticle(const G4ParticleDefinition* p)
{
  	particle = p;
  	mass = particle->GetPDGMass();
  	spin = particle->GetPDGSpin();
  	if(0.0 != spin) { spin = 0.5; }
  	tkin = 0.0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline G4double G4ScreeningMottCrossSection::GetMom2CM() const
{
	return mom2;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline G4double G4ScreeningMottCrossSection::GetMom2Lab() const
{
        return momLab2;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ScreeningMottCrossSection::GetTrec() const
{
        return Trec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ScreeningMottCrossSection::GetScreeningCoefficient() const
{
	return As;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ScreeningMottCrossSection::GetTotalCross() const
{
	return TotalCross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
