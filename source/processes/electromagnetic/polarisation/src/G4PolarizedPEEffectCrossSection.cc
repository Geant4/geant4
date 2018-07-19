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
// $Id: G4PolarizedPEEffectCrossSection.cc 96114 2016-03-16 18:51:33Z gcosmo $
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedPEEffectCrossSection
//
// Author:        Karim Laihem 
//
// Creation date: 15.03.2007
//
// Modifications:
//   19-03-07 Modified to fit in g4.8.2 framework (A.Schaelicke)
//
// Class Description:
//
#include "G4PolarizedPEEffectCrossSection.hh"
#include "G4PhysicalConstants.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedPEEffectCrossSection::G4PolarizedPEEffectCrossSection()
  {
   cout<<"G4PolarizedPEEffectCrossSection() init\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedPEEffectCrossSection::~G4PolarizedPEEffectCrossSection()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedPEEffectCrossSection::Initialize(G4double aGammaE, 
						 G4double aLept0E, 
						 G4double sinTheta,
						 const G4StokesVector & beamPol,
						 const G4StokesVector & /*p1*/,
						 G4int /*flag*/)
{
  //   cout<<"G4PolarizedPEEffectCrossSection::Initialize()\n";

// G4StokesVector PolarizedPhotoElectricEffect::Transfer_G4StokesVector(
//     G4double aGammaE,                 // Incoming Primary Gamma Energy.
//     G4ThreeVector aGammaDir,          // Incoming Primary Gamma Direction.
//     G4StokesVector beamPol,          // Incoming Primary Gamma  polarization.
//     G4double aLept0E,                 // The Lepton  e- of interest Total energy.
//     G4ThreeVector aParticl_01_Dir,    // The Lepton  e- of interest direction.
//     G4double cos_aTetha_Angle            // The lepton of interest Scattering angle.
//     )

// ***********************************************************
// ************  added by Karim   Polarization transfer to e- in PhotoelectricEffect.
// ************ 
// ***********************************************************
    G4double Gfactor   = aLept0E/electron_mass_c2+1.;
    G4double Gfactor_2 = Gfactor * Gfactor;

    G4double BETA    = sqrt(1. - 1./(Gfactor_2));

    G4double Stokes_P3  = beamPol.z()   ;                       
    
    G4double m0_c2  = electron_mass_c2; 
    G4double Lept0E = aLept0E/m0_c2+1.,   Lept0E2 = Lept0E * Lept0E ;
    G4double GammaE = aGammaE/m0_c2;


//     G4double cosTheta = cos_aTetha_Angle;
//     G4double sinTheta = sqrt(1- cos_aTetha_Angle * cos_aTetha_Angle);
    G4double cosTheta = std::sqrt(1. - sinTheta*sinTheta);

    G4double D_Lepton0 = (1./GammaE) * ((2./(GammaE*Lept0E*(1-BETA*cosTheta)))-1.);   

    G4double I_Lepton0 = 1.0+D_Lepton0;

    G4double A_Lepton0 = (Lept0E/(Lept0E+1))*(2.0/(GammaE*Lept0E) 
					      + BETA*cosTheta 
					      +(2.0/((GammaE*Lept0E2)*(1.0-BETA*cosTheta)))) / I_Lepton0 ;
    
    G4double B_Lepton0 = (Lept0E/(Lept0E+1.0)) * BETA * sinTheta * (2.0/(GammaE*Lept0E*(1-BETA*cosTheta))-1.0)/I_Lepton0;   

    G4double Stokes_S1 = (Stokes_P3 * B_Lepton0) ;
    G4double Stokes_S2 = 0.;
    G4double Stokes_S3 = (Stokes_P3 * A_Lepton0) ; 

    
    theFinalElectronPolarization.setX(Stokes_S1);
    theFinalElectronPolarization.setY(Stokes_S2);
    theFinalElectronPolarization.setZ(Stokes_S3);

    if((theFinalElectronPolarization.x()*theFinalElectronPolarization.x()
    	+ theFinalElectronPolarization.y()* theFinalElectronPolarization.y()
	+ theFinalElectronPolarization.z()* theFinalElectronPolarization.z())>1)
	
    {
	cout<<"Warning: PhotoelectricEffect Problem in pol-transfer photon to lepton:Px2 + Py2 + Pz2 > 1"<<endl;
	cout<<"Polarization transfer forced to be total and similar as incoming Photo"<<endl;
	// *KL* Surprising if it arrives (never seen it up to now)
	theFinalElectronPolarization = beamPol; // suplement de securite
// 	cout<<"PhotoEffect okay :"
//  	      <<"\t"<<(aLept0E-m0_c2)/aGammaE
// 	      <<"\t"<<aGammaE
// 	      <<"\t"<<aLept0E
// 	      <<"\t"<<cos_aTetha_Angle
// 	      <<"\t"<<beamPol
// 	      <<"\t"<<theFinalElectronPolarization
// 	      <<"\t"<<A_Lepton0
// 	      <<endl;
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4PolarizedPEEffectCrossSection::XSection(const G4StokesVector & /*pol2*/,
						   const G4StokesVector & /*pol3*/)
{
  cout<<"ERROR dummy routine G4PolarizedPEEffectCrossSection::XSection() called\n";
  return 0.;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StokesVector G4PolarizedPEEffectCrossSection::GetPol2()
{
  return theFinalElectronPolarization;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StokesVector G4PolarizedPEEffectCrossSection::GetPol3()
{
  return G4StokesVector();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
