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
// $Id: G4PolarizedMollerCrossSection.cc,v 1.1 2006-09-21 21:35:11 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

G4PolarizedMollerCrossSection::G4PolarizedMollerCrossSection()
{
}
G4PolarizedMollerCrossSection::~G4PolarizedMollerCrossSection() {}
void G4PolarizedMollerCrossSection::Initialize(
			  G4double e,
			  G4double gamma,
			  G4double phi,
		    const G4StokesVector & pol0,
                    const G4StokesVector & pol1,
			  int flag)
{
  G4double re2 = classic_electr_radius * classic_electr_radius;
  G4double gamma2=gamma*gamma;
  G4double gmo2 = (gamma - 1.)*(gamma - 1.);
  G4double pref = gamma2*re2/(gmo2*(gamma + 1.0));
  G4double sqrttwo=sqrt(2.);
  G4double w = e*(1. - e);
  G4double coeff = 0.;

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
    G4double xx = 0.;
	xx += usephi*cos(2.*phi)*((1./gamma2 - 1.) + (gamma - 1.)/(2.*gamma2*w));
	xx += (1. - 3.*gamma)/(2.*gamma2*w) - 2.*(gamma - 1.)/gamma2;
	xx *= 0.25;
    G4double yy = 0.;
	yy += -usephi*cos(2.*phi)*((1./gamma2 - 1.) + (gamma - 1.)/(2.*gamma2*w));
	yy += (1. - 3.*gamma)/(2.*gamma2*w) - 2.*(gamma - 1.)/gamma2;
	yy *= 0.25;
    G4double zz = 0.;
	zz += (1. - 2.*gamma)/(gamma*w);
	zz += (1. + 2./gamma - 3./gamma2);
	zz *= 0.25;
    phi0 += xx*pol0.x()*pol1.x() + yy*pol0.y()*pol1.y() + zz*pol0.z()*pol1.z();

    if (usephi==1.) {
    //    G4cout<<"Non-diagonal Matrix Elements"<<G4endl;
      G4double xy, yx;
	xy = -sin(2.*phi)*((1./gamma2 - 1.) + (gamma - 1.)/(2.*gamma2*w));
	xy *= 0.25;
	yx = xy;

      G4double xz, zx, zy, yz;
		coeff = sqrt((gamma + 1.)/w);
		coeff *= (2.*e - 1.)*(gamma - 1.)/(2.*sqrttwo*gamma2);
	xz = zx = -coeff*cos(phi);
	yz = zy = coeff*sin(phi);
      phi0+=yx*pol0.y()*pol1.x() + xy*pol0.x()*pol1.y();
      phi0+=zx*pol0.z()*pol1.x() + xz*pol0.x()*pol1.z();
      phi0+=zy*pol0.z()*pol1.y() + yz*pol0.y()*pol1.z();
    }
  }
  // Final state polarisarion dependence
  phi2=G4ThreeVector();
  phi3=G4ThreeVector();

  if (flag>=1) {
    G4double q = gamma*e*(e*(gamma - 1.) + 2.);
    G4double r = e*gamma*(e*(gamma - 1.) - gamma - 1.);
    G4double w2 = w*w;
    G4double sq12 = sqrt((gamma + 1.)*w);
    G4double sq32 = sqrt((gamma + 1.)*w)*(gamma + 1.)*w;
    //
    // Final Electron P1
    //

    // initial electron K1
    if (!pol0.IsZero()) {
      coeff  = 2.*e*(gamma2*gamma + (2.*w + 1.)*gamma2 + gamma*(3. - 2.*w) - 1.);
      coeff += 2.*gamma2*gamma*(w - 1.) + gamma2*(3.*w - 2.) - 6.*w*gamma + w;
      G4double xxP1K1 = coeff*cos(2.*phi);
               xxP1K1 += (gamma - 1.)*(4.*gamma*w*e - 2.*e + gamma*w + w);
               xxP1K1 /= 8.*gamma*r*w;
      G4double yyP1K1 = coeff*cos(2.*phi);
               yyP1K1 += -(gamma - 1.)*(4.*gamma*w*e - 2.*e + gamma*w + w);
	       yyP1K1 /= 8.*gamma*r*w;
      G4double zzP1K1 = (w - 1.)*gamma2 - (2.*w + 1.)*gamma + w;
	       zzP1K1 += e*((-2.*w + 1.)*gamma2 + 5.*gamma + 2.*w - 2.);
	       zzP1K1 /= 4.*r*w;
      G4double xyP1K1 = -coeff*sin(2.*phi);
	       xyP1K1 /= 8.*gamma*r*w;
      G4double yxP1K1 = -xyP1K1;
			coeff = (gamma2 - 1.)*(3.*e - 2.*w + 1.);
			coeff /= 2.*sqrttwo*r*sq12;
      G4double xzP1K1 = coeff*cos(phi);
      G4double yzP1K1 = coeff*sin(phi);
			coeff = (2.*e + 1.)*(gamma2 - 1.)*w;
			coeff /= 2.*sqrttwo*e*r*sq12;
      G4double zxP1K1 = -coeff*cos(phi);
      G4double zyP1K1 = coeff*sin(phi);
      phi2[0] += xxP1K1*pol0.x() + xyP1K1*pol0.y() + xzP1K1*pol0.z();
      phi2[1] += yxP1K1*pol0.x() + yyP1K1*pol0.y() + yzP1K1*pol0.z();
      phi2[2] += zxP1K1*pol0.x() + zyP1K1*pol0.y() + zzP1K1*pol0.z();
    }
    // initial electron K2
    if (!pol1.IsZero()) {
      coeff  = w*((3.*w + 4.)*gamma2 - 4.*w*gamma + w);
      coeff += e*(4.*(w2 + w - 1.)*gamma2 + 2.*(3. - 2.*w)*w*gamma - 2.*w);
      G4double xxP1K2 =coeff*e*cos(2.*phi);
	       xxP1K2 += (gamma - 1.)*w*((2.*gamma2 + gamma - 1.)*e*e +
					  (-2.*w*gamma2 + (3.*w - 1.)*gamma + w - 1.)*e + 
					  w*(-4.*w*gamma + gamma +1.));
	       xxP1K2 /= 8.*e*gamma*r*w2;
      G4double yyP1K2 =coeff*e*cos(2.*phi);
	       yyP1K2 += -(gamma - 1.)*w*((2.*gamma2 + gamma - 1.)*e*e + 
					   (-2.*w*gamma2 + (3.*w - 1.)*gamma + w - 1.)*e + 
					   w*(-4.*w*gamma + gamma +1.));
	       yyP1K2 /= 8.*e*gamma*r*w2;
      G4double zzP1K2 = w*(w*gamma2 + 2.*(w + 1.)*gamma - 3.*w);
	       zzP1K2 += e*((w - 2.*w2)*gamma2 + (w - 2.)*gamma + 2.*w*(w + 1.));
	       zzP1K2 /= 4.*r*w2;
      G4double xyP1K2 = -coeff*sin(2.*phi);
               xyP1K2 /= 8.*gamma*r*w2;
      G4double yxP1K2 = -xyP1K2;
			coeff = (gamma - 1.)*(gamma + 1.)*(gamma + 1.);
			coeff *= (e*(w - 2.) - 2.*w*(w - 1.));
			coeff /= 2.*sqrttwo*r*sq32;
      G4double xzP1K2 = coeff*cos(phi);
      G4double yzP1K2 = coeff*sin(phi);
			coeff = (gamma - 1.)*(gamma + 1.)*(gamma + 1.);
			coeff *= (e*(w - 2.) + 2.*w*(w + 1.));
			coeff /= 2.*sqrttwo*r*sq32;
      G4double zxP1K2 = -coeff*cos(phi);
      G4double zyP1K2 = coeff*sin(phi);
      phi2[0] += xxP1K2*pol1.x() + xyP1K2*pol1.y() + xzP1K2*pol1.z();
      phi2[1] += yxP1K2*pol1.x() + yyP1K2*pol1.y() + yzP1K2*pol1.z();
      phi2[2] += zxP1K2*pol1.x() + zyP1K2*pol1.y() + zzP1K2*pol1.z();
    }
    //
    // Final Electron P2
    //

    // initial electron K1
    if (!pol0.IsZero()) {
               coeff = (w*(2.*gamma2*gamma*(w - 1.) + (3.*w - 2.)*gamma2 - 6.*w*gamma + w) +
			e*((2. - 4.*w)*gamma2*gamma + (-4.*w2 - 5.*w + 2.)*gamma2 + 4.*w2*gamma + w));
      G4double xxP2K1  = coeff*e*cos(2.*phi);
               xxP2K1 += -(gamma - 1.)*w*(e - w)*(e*gamma + e + 4.*gamma*w - 2.);
	       xxP2K1 /= 8.*e*gamma*q*w2;
      G4double yyP2K1  = coeff*e*cos(2.*phi);
	       yyP2K1 += (gamma - 1.)*w*(e - w)*(e*gamma + e + 4.*gamma*w - 2.);
	       yyP2K1 /= 8.*e*gamma*q*w2;
      G4double zzP2K1 = w*((w - 1.)*gamma2 - (2.*w + 1.)*gamma + w);
	       zzP2K1 += e*((2.*w2 - 2.*w + 1.)*gamma2 - 3.*w*gamma + gamma -2.*w2 +w);
	       zzP2K1 /= 4.*q*w2;
      G4double xyP2K1 = -coeff*sin(2.*phi);
	       xyP2K1 /= 8.*gamma*q*w2;
      G4double yxP2K1 = -xyP2K1;
               coeff = (gamma - 1.)*(gamma + 1.)*(gamma + 1.)*(e*(w + 1.) + w*(2.*w - 1.));
	       coeff /=  2.*sqrttwo*q*sq32;
      G4double xzP2K1 = coeff*cos(phi);
      G4double yzP2K1 = coeff*sin(phi);
               coeff = (gamma - 1.)*(gamma + 1.)*(gamma + 1.)*(e*(w + 1.) - w*(2.*w + 1.));
	       coeff /=  2.*sqrttwo*q*sq32;
      G4double zxP2K1 = -coeff*cos(phi);
      G4double zyP2K1 = coeff*sin(phi);
      phi3[0] += xxP2K1*pol0.x() + xyP2K1*pol0.y() + xzP2K1*pol0.z();
      phi3[1] += yxP2K1*pol0.x() + yyP2K1*pol0.y() + yzP2K1*pol0.z();
      phi3[2] += zxP2K1*pol0.x() + zyP2K1*pol0.y() + zzP2K1*pol0.z();
    }
	// initial electron K2
    if (!pol1.IsZero()) {
      coeff  = e*(-(4.*w + 7.)*gamma2 + (4.*w - 2.)*gamma + 1.);
      coeff += (3.*w + 4.)*gamma2 - (4.*gamma - 1.)*w;
      G4double xxP2K2  = coeff*cos(2.*phi);
	       xxP2K2 += -(gamma - 1.)*((2.*gamma2 + gamma -1.)*w + e*(gamma*(4.*w - 1.)- 1.));
	       xxP2K2 /= 8.*gamma*q*w;
      G4double yyP2K2  = coeff*cos(2.*phi);
	       yyP2K2 += (gamma - 1.)*((2.*gamma2 + gamma -1.)*w + e*(gamma*(4.*w - 1.)- 1.));
	       yyP2K2 /= 8.*gamma*q*w;
      G4double zzP2K2 = w*gamma2 + 2.*(w + 1.)*gamma - 3.*w + e*(2.*(w - 1.)*gamma2 
								 - 3.*gamma - 2.*w + 1.);
	       zzP2K2 /= 4.*q*w;
      G4double xyP2K2 = -coeff*sin(2.*phi);
	       xyP2K2 /= 8.*gamma*q*w;
      G4double yxP2K2 = -xyP2K2;
	                coeff = (gamma2 - 1.)*(e - 2.*w + 2.);
			coeff /= 2.*sqrttwo*q*sq12;
      G4double xzP2K2 = -coeff*cos(phi);
      G4double yzP2K2 = -coeff*sin(phi);
			coeff = (gamma2 - 1.)*(-3.*e + 2.*w + 2.);
			coeff /= 2.*sqrttwo*q*sq12;
      G4double zxP2K2 = coeff*cos(phi);
      G4double zyP2K2 = -coeff*sin(phi);
      phi3[0] += xxP2K2*pol1.x() + xyP2K2*pol1.y() + xzP2K2*pol1.z();
      phi3[1] += yxP2K2*pol1.x() + yyP2K2*pol1.y() + yzP2K2*pol1.z();
      phi3[2] += zxP2K2*pol1.x() + zyP2K2*pol1.y() + zzP2K2*pol1.z();
    }
  }
  phi0 *= pref;
  phi2 *= pref;
  phi3 *= pref;
}

G4double G4PolarizedMollerCrossSection::XSection(const G4StokesVector & pol2,const G4StokesVector & pol3)
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
  const G4double xmin, const G4double xmax, const G4double gamma,
  const G4StokesVector & pol0,const G4StokesVector & pol1)
{
  G4double xs=0.;

  G4double x=xmin;

  if (xmax != 1./2.) G4cout<<" warning xmax expected to be 1/2 but is "<<xmax<< G4endl;

  //  re -> electron radius^2;
  G4double re2 = classic_electr_radius * classic_electr_radius;
  G4double gamma2=gamma*gamma;
  G4double gmo2 = (gamma - 1.)*(gamma - 1.);
  G4double logMEM = log(1./x - 1.);
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
