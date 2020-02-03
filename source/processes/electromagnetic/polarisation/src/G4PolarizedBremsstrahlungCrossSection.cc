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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedBremsstrahlungCrossSection
//
// Author:        Andreas Schaelicke on the base of Karim Laihems code
//
// Creation date: 16.08.2006
//

#include "G4PolarizedBremsstrahlungCrossSection.hh"
#include "G4PhysicalConstants.hh"

G4bool G4PolarizedBremsstrahlungCrossSection::scrnInitialized=false;
G4double G4PolarizedBremsstrahlungCrossSection::SCRN [3][20];  
// screening function lookup table;


void G4PolarizedBremsstrahlungCrossSection::InitializeMe()
{
  if (!scrnInitialized) {
    SCRN [1][1]=  0.5   ; SCRN [2][1] = 0.0145;
    SCRN [1][2]=  1.0   ; SCRN [2][2] = 0.0490;
    SCRN [1][3]=  2.0   ; SCRN [2][3] = 0.1400;
    SCRN [1][4]=  4.0   ; SCRN [2][4] = 0.3312;
    SCRN [1][5]=  8.0   ; SCRN [2][5] = 0.6758;
    SCRN [1][6]=  15.0  ; SCRN [2][6] = 1.126;
    SCRN [1][7]=  20.0  ; SCRN [2][7] = 1.367;
    SCRN [1][8]=  25.0  ; SCRN [2][8] = 1.564;
    SCRN [1][9]=  30.0  ; SCRN [2][9] = 1.731;
    SCRN [1][10]= 35.0  ; SCRN [2][10]= 1.875;
    SCRN [1][11]= 40.0  ; SCRN [2][11]= 2.001;
    SCRN [1][12]= 45.0  ; SCRN [2][12]= 2.114;
    SCRN [1][13]= 50.0  ; SCRN [2][13]= 2.216;
    SCRN [1][14]= 60.0  ; SCRN [2][14]= 2.393;
    SCRN [1][15]= 70.0  ; SCRN [2][15]= 2.545;
    SCRN [1][16]= 80.0  ; SCRN [2][16]= 2.676;
    SCRN [1][17]= 90.0  ; SCRN [2][17]= 2.793;
    SCRN [1][18]= 100.0 ; SCRN [2][18]= 2.897;
    SCRN [1][19]= 120.0 ; SCRN [2][19]= 3.078;   

    scrnInitialized=true;
  }
}

G4PolarizedBremsstrahlungCrossSection::G4PolarizedBremsstrahlungCrossSection()
{
  InitializeMe();
}


void G4PolarizedBremsstrahlungCrossSection::Initialize(
     G4double aLept0E, G4double aGammaE, G4double sintheta,
     const G4StokesVector & beamPol,
     const G4StokesVector & /*p1*/,
     G4int /*flag*/)
{
//   G4cout<<"G4PolarizedBremsstrahlungCrossSection::Initialize \n"
// 	<<"lepE = "<<aLept0E
// 	<<"gamE = "<<aGammaE
// 	<<"sint = "<<sintheta<<"\n"
// 	<<"beamPol="<<beamPol<<"\n";

  G4double aLept1E = aLept0E - aGammaE;

  G4double Stokes_S1  = beamPol.x()   ; 
  G4double Stokes_S2  = beamPol.y()   ; 
  G4double Stokes_S3  = beamPol.z()   ; 
  // **************************************************************************  
    
  G4double m0_c2  = electron_mass_c2; 
  G4double Lept0E = aLept0E/m0_c2+1.,   Lept0E2 = Lept0E * Lept0E ;
  G4double GammaE = aGammaE/m0_c2,      GammaE2 = GammaE * GammaE ;
  G4double Lept1E = aLept1E/m0_c2+1.,   Lept1E2 = Lept1E * Lept1E ;
    

  //  const G4Element* theSelectedElement = theModel->SelectedAtom();

  // *******  Gamma Transvers Momentum
    
  G4double TMom = std::sqrt(Lept0E2 -1.)* sintheta;
  G4double u    = TMom       , u2 =u * u ;
  G4double Xsi  = 1./(1.+u2)                      , Xsi2 = Xsi * Xsi  ; 

  //  G4double theZ  = theSelectedElement->GetZ();
    
  //  G4double fCoul = theSelectedElement->GetfCoulomb();
  G4double delta = 12. * std::pow(theZ, 1./3.) * 
    Lept0E * Lept1E * Xsi / (121. * GammaE); 
  G4double GG=0.;

  if(delta < 0.5) {
    GG = std::log(2.* Lept0E * Lept1E / GammaE) - 2. - fCoul; 
  }
  else if ( delta < 120) {
    for (G4int j=2; j<=19; j++)  {
      if(SCRN[1][j] >= delta)    {
	GG =std::log(2 * Lept0E * Lept1E / GammaE) - 2 - fCoul
	  -(SCRN[2][j-1]+(delta-SCRN[1][j-1])*(SCRN[2][j]-SCRN[2][j-1])
	    /(SCRN[1][j]-SCRN[1][j-1]));
	break;
      }
    }
  }
  else  {
    G4double alpha_sc  = (111. * std::pow(theZ, -1./3.)) / Xsi;
    GG = std::log(alpha_sc)- 2. - fCoul;
  }

  if(GG<-1.) GG=-1.;     // *KL* do we need this ?!

  G4double I_Lept   = (Lept0E2 + Lept1E2) * (3.+2.*GG) - 2 * Lept0E * Lept1E * (1. + 4. * u2 * Xsi2 * GG);
  G4double F_Lept   = Lept1E * 4. * GammaE *  u * Xsi * (1. - 2 * Xsi) * GG / I_Lept;
  G4double E_Lept   = Lept0E * 4. * GammaE *  u * Xsi * (2. * Xsi - 1.) * GG / I_Lept; 
  G4double M_Lept   = 4. * Lept0E * Lept1E * (1. + GG - 2. * Xsi2 * u2 * GG) / I_Lept ;
  G4double P_Lept   = GammaE2 * (1. + 8. * GG * (Xsi - 0.5)*(Xsi - 0.5)) / I_Lept ;
    
  G4double Stokes_SS1 = M_Lept * Stokes_S1 + E_Lept * Stokes_S3;
  G4double Stokes_SS2 = M_Lept * Stokes_S2 ;
  G4double Stokes_SS3 = (M_Lept + P_Lept) * Stokes_S3 + F_Lept * Stokes_S1; 

  theFinalLeptonPolarization.setX(Stokes_SS1);
  theFinalLeptonPolarization.setY(Stokes_SS2);
  theFinalLeptonPolarization.setZ(Stokes_SS3);

  if(theFinalLeptonPolarization.mag2()>1) { 
    G4cout<<" WARNING in pol-brem theFinalLeptonPolarization \n";
    G4cout
      <<"\t"<<theFinalLeptonPolarization
      <<"\t GG\t"<<GG
      <<"\t delta\t"<<delta
      <<G4endl;
    theFinalLeptonPolarization.setX(0);
    theFinalLeptonPolarization.setY(0);
    theFinalLeptonPolarization.setZ(Stokes_SS3);
    if(Stokes_SS3>1) theFinalLeptonPolarization.setZ(1);
  }


  G4double I_Gamma   = (Lept0E2 + Lept1E2)*(3.+2.*GG) - 2. * Lept0E * Lept1E * (1. + 4. * u2 * Xsi2 * GG);
  G4double D_Gamma   = 8. * Lept0E * Lept1E * u2 * Xsi2 * GG / I_Gamma;
  G4double L_Gamma   = GammaE * ((Lept0E + Lept1E) * (3. + 2. * GG) 
				 - 2. * Lept1E * (1. + 4. * u2 * Xsi2 * GG))/I_Gamma;   
  G4double T_Gamma   = 4. * GammaE * Lept1E * Xsi * u * (2. * Xsi - 1.) * GG / I_Gamma ;
    
  G4double Stokes_P1 = D_Gamma ;
  G4double Stokes_P2 = 0.;
  G4double Stokes_P3 = (Stokes_S3*L_Gamma + Stokes_S1*T_Gamma) ;

  theFinalGammaPolarization.SetPhoton();

  theFinalGammaPolarization.setX(Stokes_P1);
  theFinalGammaPolarization.setY(Stokes_P2);
  theFinalGammaPolarization.setZ(Stokes_P3);
    
  if(theFinalGammaPolarization.mag2()>1) {
    G4cout<<" WARNING in pol-brem theFinalGammaPolarization \n";
    G4cout
      <<"\t"<<theFinalGammaPolarization
      <<"\t GG\t"<<GG
      <<"\t delta\t"<<delta
      <<G4endl;
  }
}

G4double G4PolarizedBremsstrahlungCrossSection::XSection(const G4StokesVector & /*pol2*/,
							 const G4StokesVector & /*pol3*/)
{
  G4cout<<"ERROR dummy routine G4PolarizedBremsstrahlungCrossSection::XSection called \n";
  return 0.;
}

  // return expected mean polarisation
G4StokesVector G4PolarizedBremsstrahlungCrossSection::GetPol2()
{
  // electron/positron
  return theFinalLeptonPolarization;
}
G4StokesVector G4PolarizedBremsstrahlungCrossSection::GetPol3()
{
  // photon
  return theFinalGammaPolarization;;
}


