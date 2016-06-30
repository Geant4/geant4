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
// $Id: G4PolarizedComptonCrossSection.cc 96114 2016-03-16 18:51:33Z gcosmo $
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedComptonCrossSection
//
// Author:        Andreas Schaelicke
//
// Creation date: 15.05.2005
//
// Modifications:
//
// Class Description:
//   determine the  polarization of the final state 
//   in a Compton scattering process employing the differential 
//   cross section by F.W.Lipps & H.A.Tolhoek
//   ( Physica 20 (1954) 395 )
//   recalculated by P.Starovoitov
//
#include "G4PolarizedComptonCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedComptonCrossSection::G4PolarizedComptonCrossSection() 
  : gammaPol2(false), electronPol3(false)
{
  SetYmin(0.);

  //  G4cout<<"G4PolarizedComptonCrossSection() init\n";

  re2 = classic_electr_radius * classic_electr_radius * sqr(4*pi/hbarc);
  //  G4double unit_conversion = hbarc_squared ;
  //  G4cout<<" (keV)^2* m^2 ="<<unit_conversion<<"\n";
  phi0 = 0.; polXS = 0.; unpXS = 0.;
  phi2 = G4ThreeVector(0., 0., 0.);
  phi3 = G4ThreeVector(0., 0., 0.);
  polxx = polyy = polzz = polxz = polzx = polyz = polzy = polxy = polyx = 0.;
  diffXSFactor = 1.;
  totalXSFactor = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedComptonCrossSection::~G4PolarizedComptonCrossSection()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedComptonCrossSection::Initialize(G4double eps, G4double X, G4double , // phi
						const G4StokesVector & pol0,
						const G4StokesVector & pol1,
						G4int flag)
{
  G4double cosT = 1. - (1./eps - 1.)/X;
  if(cosT > 1.+1.e-8)  cosT = 1.;
  if(cosT < -1.-1.e-8) cosT = -1.;
  G4double cosT2 = cosT*cosT;
  G4double cosT3 = cosT2*cosT;
  G4double sinT2 = 1. - cosT2;
  if(sinT2 > 1. + 1.e-8) sinT2 = 1.;
  if(sinT2 < 0.)     sinT2 = 0.;
  G4double sinT = std::sqrt(sinT2);
  G4double cos2T = 2.*cosT2 - 1.;
  G4double sin2T = 2.*sinT*cosT;
  G4double eps2 = sqr(eps);
  DefineCoefficients(pol0,pol1);
  diffXSFactor = re2/(4.*X);

  // unpolarized Cross Section
  unpXS = (eps2 + 1. - eps*sinT2)/(2.*eps);
  // initial polarization dependence
  polXS = -sinT2*pol0.x() + (1. - eps)*sinT*polzx + ((eps2 - 1.)/eps)*cosT*polzz;
  polXS *= 0.5;

  phi0 = unpXS + polXS;

  if (flag == 2 ){
  // polarization of outgoing photon
  G4double PHI21 = -sinT2 + 0.5*(cos2T + 3.)*pol0.x() - ((1. - eps)/eps)*sinT*polzx;
  PHI21 *= 0.5;
  G4double PHI22 = cosT*pol0.y() + ((1. - eps)/(2.*eps))*sinT*polzy;
  G4double PHI23 = ((eps2 + 1.)/eps)*cosT*pol0.z() - ((1. - eps)/eps)*(eps*cosT2 + 1.)*pol1.z();
  PHI23 += 0.5*(1. - eps)*sin2T*pol1.x();
  PHI23 += (eps - 1.)*(-sinT2*polxz + sinT*polyy - 0.5*sin2T*polxx);
  PHI23 *= 0.5;
  phi2 = G4ThreeVector(PHI21, PHI22, PHI23);

  // polarization of outgoing electron
  G4double PHI32 = -sinT2*polxy + ((1. - eps)/eps)*sinT*polyz + 0.5*(cos2T + 3.)*pol1.y();
  PHI32 *= 0.5;

  G4double PHI31 = 0., PHI31add = 0., PHI33 = 0., PHI33add = 0.;

  if ((1. - eps) > 1.e-12){
    G4double helpVar = std::sqrt(eps2 - 2.*cosT*eps + 1.);

    PHI31 = (1. - eps)*(1. + cosT)*sinT*pol0.z();
    PHI31 += (-eps*cosT3 + eps*cosT2 + (eps - 2.)*cosT + eps)*pol1.x();
    PHI31 += -(eps*cosT2 - eps*cosT + cosT + 1.)*sinT*pol1.z();
    PHI31 /= 2.*helpVar;

    PHI31add = -eps*sqr(1. - cosT)*(1. + cosT)*polxx;
    PHI31add += (1. - eps)*sinT2*polyy;
    PHI31add += -(-eps2 + cosT*(cosT*eps - eps + 1.)*eps + eps - 1.)*sinT*polxz/eps;
    PHI31add /= 2.*helpVar;
    
    PHI33 = ((1. - eps)/eps)*(-eps*cosT2 + eps*(eps + 1.)*cosT - 1.)*pol0.z();
    PHI33 += -(eps*cosT2 + (1. - eps)*eps*cosT + 1.)*sinT*pol1.x();
    PHI33 += -(-eps2*cosT3 + eps*(eps2 - eps + 1.)*cosT2 - cosT + eps2)*pol1.z()/eps;
    PHI33 /= -2.*helpVar;

    PHI33add = (eps*(eps - cosT - 1.)*cosT + 1.)*sinT*polxx;
    PHI33add += -(-eps2 + cosT*eps + eps - 1.)*sinT2*polxz;
    PHI33add += (eps - 1.)*(cosT - eps)*sinT*polyy;
    PHI33add /= -2.*helpVar;
  }else{
     PHI31 = -pol1.z() - (X - 1.)*std::sqrt(1. - eps)*pol1.x()/std::sqrt(2.*X);
     PHI31add = -(-X*X*pol1.z() - 2.*X*(2.*pol0.z() - pol1.z()) - (4.*pol0.x() + 5.)*pol1.z())*(1. - eps)/(4.*X);
     
     PHI33 = pol1.x() - (X - 1.)*std::sqrt(1. - eps)*pol1.z()/std::sqrt(2.*X);
     PHI33add = -(X*X - 2.*X + 4.*pol0.x() + 5.)*(1. - eps)*pol1.x()/(4.*X);
  }
  phi3 = G4ThreeVector(PHI31 + PHI31add, PHI32, PHI33 + PHI33add);
    
  }
  unpXS *= diffXSFactor;
  polXS *= diffXSFactor;
  phi0 *= diffXSFactor;
  phi2 *= diffXSFactor;
  phi3 *= diffXSFactor;
  
}

G4double G4PolarizedComptonCrossSection::XSection(const G4StokesVector & pol2,const G4StokesVector & pol3)
{
  gammaPol2    = !(pol2==G4StokesVector::ZERO);
  electronPol3 = !(pol3==G4StokesVector::ZERO);

  G4double phi = 0.;
  // polarization independent part
  phi += phi0;


  if (gammaPol2) { 
    // part depending on the polarization of the final photon  
    phi += phi2*pol2;
  }

  if (electronPol3) {
    // part depending on the polarization of the final electron  
    phi += phi3*pol3;
  }

  // return cross section.
  return phi;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4PolarizedComptonCrossSection::TotalXSection(G4double /*xmin*/, G4double /*xmax*/, G4double k0,
						   const G4StokesVector & pol0,
						   const G4StokesVector & pol1)
{
  
  //  G4double k0 = gammaEnergy / electron_mass_c2 ;
  G4double k1 = 1. + 2.*k0 ;

//   // pi*re^2
//   G4double re=2.81794e-15; //m
//   G4double barn=1.e-28; //m^2
  G4double Z=theZ;

  G4double unit = Z*pi*classic_electr_radius  
    * classic_electr_radius ; // *1./barn;

  G4double pre = unit/(sqr(k0)*sqr(1.+2.*k0));

  G4double xs_0 = ((k0 - 2.)*k0  -2.)*sqr(k1)*std::log(k1) + 2.*k0*(k0*(k0 + 1.)*(k0 + 8.) + 2.);		
  G4double xs_pol = (k0 + 1.)*sqr(k1)*std::log(k1) - 2.*k0*(5.*sqr(k0) + 4.*k0 + 1.);

  return pre*(xs_0/k0 + pol0.p3()*pol1.z()*xs_pol);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4StokesVector G4PolarizedComptonCrossSection::GetPol2()
{
  // Note, mean polarization can not contain correlation 
  // effects.
  return 1./phi0 * phi2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StokesVector G4PolarizedComptonCrossSection::GetPol3()
{
  // Note, mean polarization can not contain correlation 
  // effects.
  return 1./phi0 * phi3;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedComptonCrossSection::DefineCoefficients(const G4StokesVector & pol0,
							const G4StokesVector & pol1)
{
  polxx=pol0.x()*pol1.x();
  polyy=pol0.y()*pol1.y();
  polzz=pol0.z()*pol1.z();

  polxz=pol0.x()*pol1.z();
  polzx=pol0.z()*pol1.x();

  polyz=pol0.y()*pol1.z();
  polzy=pol0.z()*pol1.y();

  polxy=pol0.x()*pol1.y();
  polyx=pol0.y()*pol1.x();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
