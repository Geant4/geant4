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
// $Id: G4PolarizedMollerCrossSection.cc 96114 2016-03-16 18:51:33Z gcosmo $
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedMollerCrossSection
//
// Author:        Andreas Schaelicke
//
// Creation date: 12.01.2006
//
// Modifications:
//   16-01-06 included cross section as calculated by P.Starovoitov
//
// Class Description:
//   * calculates the differential cross section
//     incomming electron K1(along positive z direction) scatters at an electron K2 at rest
//   * phi denotes the angle between the scattering plane (defined by the
//     outgoing electron) and X-axis
//   * all stokes vectors refer to spins in the Global System (X,Y,Z)
//

#include "G4PolarizedMollerCrossSection.hh"
#include "G4PhysicalConstants.hh"

G4PolarizedMollerCrossSection::G4PolarizedMollerCrossSection() :
  phi0(0.)
{
  SetXmax(.5);
}
G4PolarizedMollerCrossSection::~G4PolarizedMollerCrossSection() {}
void G4PolarizedMollerCrossSection::Initialize(
					       G4double e,
					       G4double gamma,
					       G4double /*phi*/,
					       const G4StokesVector & pol0,
					       const G4StokesVector & pol1,
					       G4int flag)
{
  G4double re2 = classic_electr_radius * classic_electr_radius;
  G4double gamma2=gamma*gamma;
  G4double gmo = (gamma - 1.);
  G4double gmo2 = (gamma - 1.)*(gamma - 1.);
  G4double gpo = (gamma + 1.);
  G4double pref = gamma2*re2/(gmo2*(gamma + 1.0));
  G4double sqrttwo=std::sqrt(2.);
  G4double f = (-1. + e);
  G4double e2 = e*e;
  G4double f2 = f*f;
  //  G4double w = e*(1. - e);

  G4bool polarized=(!pol0.IsZero())||(!pol1.IsZero());

  if (flag==0) polarized=false;
  // Unpolarised part of XS
  phi0 = 0.;
  phi0+= gmo2/gamma2;
  phi0+= ((1. - 2.*gamma)/gamma2)*(1./e + 1./(1.-e));
  phi0+= 1./(e*e) + 1./((1. - e)*(1. - e));
  phi0*=0.25;
  // Initial state polarisarion dependence
  if (polarized) {
    G4double usephi=1.;
    if (flag<=1) usephi=0.;
    //    G4cout<<"Polarized differential moller cross section"<<G4endl;
    //    G4cout<<"Initial state polarisation contributions"<<G4endl;
    //    G4cout<<"Diagonal Matrix Elements"<<G4endl;
    G4double xx = (gamma - f*e*gmo*(3. + gamma))/(4.*f*e*gamma2);
    G4double yy = (-1. + f*e*gmo2 + 2.*gamma)/(4.*f*e*gamma2);
    G4double zz = (-(e*gmo*(3. + gamma)) + e2*gmo*(3. + gamma) + 
		   gamma*(-1. + 2.*gamma))/(4.*f*e*gamma2);

    phi0 += xx*pol0.x()*pol1.x() + yy*pol0.y()*pol1.y() + zz*pol0.z()*pol1.z();

    if (usephi==1.) {
      //    G4cout<<"Non-diagonal Matrix Elements"<<G4endl;
      G4double xy = 0.;
      G4double xz = -((-1. + 2.*e)*gmo)/(2.*sqrttwo*gamma2*
				       std::sqrt(-((f*e)/gpo)));
      G4double yx = 0.;
      G4double yz = 0.;
      G4double zx = -((-1. + 2.*e)*gmo)/(2.*sqrttwo*gamma2*
				       std::sqrt(-((f*e)/gpo)));
      G4double zy = 0.;
      phi0+=yx*pol0.y()*pol1.x() + xy*pol0.x()*pol1.y();
      phi0+=zx*pol0.z()*pol1.x() + xz*pol0.x()*pol1.z();
      phi0+=zy*pol0.z()*pol1.y() + yz*pol0.y()*pol1.z();
    }
  }
  // Final state polarisarion dependence
  phi2=G4ThreeVector();
  phi3=G4ThreeVector();

  if (flag>=1) {
    //
    // Final Electron P1
    //

    // initial electron K1
    if (!pol0.IsZero()) {
      G4double xxP1K1 = (std::sqrt(gpo/(1. + e2*gmo + gamma - 2.*e*gamma))*
			 (gamma - e*gpo))/(4.*e2*gamma);
      G4double xyP1K1 = 0.;
      G4double xzP1K1 = (-1. + 2.*e*gamma)/(2.*sqrttwo*f*gamma*
					  std::sqrt(e*e2*(1. + e + gamma - e*gamma)));
      G4double yxP1K1 = 0.;
      G4double yyP1K1 = (-gamma2 + e*(-1. + gamma*(2. + gamma)))/(4.*f*e2*gamma2);
      G4double yzP1K1 = 0.;
      G4double zxP1K1 = (1. + 2.*e2*gmo - 2.*e*gamma)/(2.*sqrttwo*f*e*gamma*
						    std::sqrt(e*(1. + e + gamma - e*gamma)));
      G4double zyP1K1 = 0.;
      G4double zzP1K1 = (-gamma + e*(1. - 2.*e*gmo + gamma))/(4.*f*e2*gamma*
							    std::sqrt(1. - (2.*e)/(f*gpo)));
      phi2[0] += xxP1K1*pol0.x() + xyP1K1*pol0.y() + xzP1K1*pol0.z();
      phi2[1] += yxP1K1*pol0.x() + yyP1K1*pol0.y() + yzP1K1*pol0.z();
      phi2[2] += zxP1K1*pol0.x() + zyP1K1*pol0.y() + zzP1K1*pol0.z();
    }
    // initial electron K2
    if (!pol1.IsZero()) {
      G4double xxP1K2 = ((1. + e*(-3. + gamma))*std::sqrt(gpo/(1. + e2*gmo + gamma - 
							     2.*e*gamma)))/(4.*f*e*gamma);
      G4double xyP1K2 = 0.;
      G4double xzP1K2 = (-2. + 2.*e + gamma)/(2.*sqrttwo*f2*gamma*
					    std::sqrt(e*(1. + e + gamma - e*gamma)));
      G4double yxP1K2 = 0.;
      G4double yyP1K2 = (1. - 2.*gamma + e*(-1. + gamma*(2. + gamma)))/(4.*f2*e*gamma2);
      G4double yzP1K2 = 0.;
      G4double zxP1K2 = (2.*e*(1. + e*gmo - 2.*gamma) + gamma)/(2.*sqrttwo*f2*gamma*
							     std::sqrt(e*(1. + e + gamma - e*gamma)));
      G4double zyP1K2 = 0.;
      G4double zzP1K2 = (1. - 2.*gamma + e*(-1. - 2.*e*gmo + 3.*gamma))/
	(4.*f2*e*gamma*std::sqrt(1. - (2.*e)/(f*gpo)));
      phi2[0] += xxP1K2*pol1.x() + xyP1K2*pol1.y() + xzP1K2*pol1.z();
      phi2[1] += yxP1K2*pol1.x() + yyP1K2*pol1.y() + yzP1K2*pol1.z();
      phi2[2] += zxP1K2*pol1.x() + zyP1K2*pol1.y() + zzP1K2*pol1.z();
    }
    //
    // Final Electron P2
    //

    // initial electron K1
    if (!pol0.IsZero()) {


      G4double xxP2K1 = (-1. + e + e*gamma)/(4.*f2*gamma*
					    std::sqrt((e*(2. + e*gmo))/gpo));
      G4double xyP2K1 = 0.;
      G4double xzP2K1 = -((1. + 2.*f*gamma)*std::sqrt(f/(-2. + e - e*gamma)))/
	(2.*sqrttwo*f2*e*gamma);
      G4double yxP2K1 = 0.;
      G4double yyP2K1 = (1. - 2.*gamma + e*(-1. + gamma*(2. + gamma)))/(4.*f2*e*gamma2);
      G4double yzP2K1 = 0.;
      G4double zxP2K1 = (1. + 2.*e*(-2. + e + gamma - e*gamma))/(2.*sqrttwo*f*e*
							      std::sqrt(-(f*(2. + e*gmo)))*gamma);
      G4double zyP2K1 = 0.;
      G4double zzP2K1 = (std::sqrt((e*gpo)/(2. + e*gmo))*
			 (-3. + e*(5. + 2.*e*gmo - 3.*gamma) + 2.*gamma))/(4.*f2*e*gamma);

      phi3[0] += xxP2K1*pol0.x() + xyP2K1*pol0.y() + xzP2K1*pol0.z();
      phi3[1] += yxP2K1*pol0.x() + yyP2K1*pol0.y() + yzP2K1*pol0.z();
      phi3[2] += zxP2K1*pol0.x() + zyP2K1*pol0.y() + zzP2K1*pol0.z();
    }
    // initial electron K2
    if (!pol1.IsZero()) {

      G4double xxP2K2 = (-2. - e*(-3. + gamma) + gamma)/
	(4.*f*e*gamma* std::sqrt((e*(2. + e*gmo))/gpo));
      G4double xyP2K2 = 0.;
      G4double xzP2K2 = ((-2.*e + gamma)*std::sqrt(f/(-2. + e - e*gamma)))/
	(2.*sqrttwo*f*e2*gamma);
      G4double yxP2K2 = 0.;
      G4double yyP2K2 = (-gamma2 + e*(-1. + gamma*(2. + gamma)))/(4.*f*e2*gamma2);
      G4double yzP2K2 = 0.;
      G4double zxP2K2 = (gamma + 2.*e*(-1. + e - e*gamma))/
	(2.*sqrttwo*e2* std::sqrt(-(f*(2. + e*gmo)))*gamma);
      G4double zyP2K2 = 0.;
      G4double zzP2K2 = (std::sqrt((e*gpo)/(2. + e*gmo))*
			 (-2. + e*(3. + 2.*e*gmo - gamma) + gamma))/(4.*f*e2*gamma);
      phi3[0] += xxP2K2*pol1.x() + xyP2K2*pol1.y() + xzP2K2*pol1.z();
      phi3[1] += yxP2K2*pol1.x() + yyP2K2*pol1.y() + yzP2K2*pol1.z();
      phi3[2] += zxP2K2*pol1.x() + zyP2K2*pol1.y() + zzP2K2*pol1.z();
    }
  }
  phi0 *= pref;
  phi2 *= pref;
  phi3 *= pref;
}

G4double G4PolarizedMollerCrossSection::XSection(const G4StokesVector & pol2,
						 const G4StokesVector & pol3)
{
  G4double xs=0.;
  xs+=phi0;

  G4bool polarized=(!pol2.IsZero())||(!pol3.IsZero());
  if (polarized) {
    xs+=phi2*pol2 + phi3*pol3;
  }
  return xs;
}

G4double G4PolarizedMollerCrossSection::TotalXSection(
  G4double xmin, G4double xmax, G4double gamma,
  const G4StokesVector & pol0,const G4StokesVector & pol1)
{
  G4double xs=0.;

  G4double x=xmin;

  if (xmax != 1./2.) G4cout<<" warning xmax expected to be 1/2 but is "<<xmax<< G4endl;

  //  re -> electron radius^2;
  G4double re2 = classic_electr_radius * classic_electr_radius;
  G4double gamma2=gamma*gamma;
  G4double gmo2 = (gamma - 1.)*(gamma - 1.);
  G4double logMEM = std::log(1./x - 1.);
  G4double pref = twopi*gamma2*re2/(gmo2*(gamma + 1.0));
  // unpolarise XS
  G4double sigma0 = 0.;
  	sigma0 += (gmo2/gamma2)*(0.5 - x);
	sigma0 += ((1. - 2.*gamma)/gamma2)*logMEM;
	sigma0 += 1./x - 1./(1. - x);
  //    longitudinal part
  G4double sigma2=0.;
  	sigma2 += ((gamma2 + 2.*gamma - 3.)/gamma2)*(0.5 - x);
	sigma2 += (1./gamma - 2.)*logMEM;
  //    transverse part
  G4double sigma3=0.;
  	sigma3 += (2.*(1. - gamma)/gamma2)*(0.5 - x);
	sigma3 += (1. - 3.*gamma)/(2.*gamma2)*logMEM;
  // total cross section
  xs+=pref*(sigma0 + sigma2*pol0.z()*pol1.z() + sigma3*(pol0.x()*pol1.x()+pol0.y()*pol1.y()));

  return xs;
}


G4StokesVector G4PolarizedMollerCrossSection::GetPol2()
{
  // Note, mean polarization can not contain correlation
  // effects.
  return 1./phi0 * phi2;
}
G4StokesVector G4PolarizedMollerCrossSection::GetPol3()
{
  // Note, mean polarization can not contain correlation
  // effects.
  return 1./phi0 * phi3;
}
