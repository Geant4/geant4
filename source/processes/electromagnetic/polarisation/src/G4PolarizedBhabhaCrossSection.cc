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
// $Id: G4PolarizedBhabhaCrossSection.cc,v 1.2 2006-09-26 09:08:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PolarizedBhabhaCrossSection
//
// Author:        Andreas Schaelicke
//
// Creation date: 12.01.2006
//
// Modifications:
//   16-01-06 included cross section as calculated by P.Starovoitov
//   24-08-06 bugfix in total cross section (A. Schaelicke)
//
// Class Description:
//   * calculates the differential cross section
//     incomming positron Kpl(along positive z direction) scatters at an electron Kmn at rest
//   * phi denotes the angle between the scattering plane (defined by the
//     outgoing electron) and X-axis
//   * all stokes vectors refer to spins in the Global System (X,Y,Z)
//

#include "G4PolarizedBhabhaCrossSection.hh"

G4PolarizedBhabhaCrossSection::G4PolarizedBhabhaCrossSection()
{
}
G4PolarizedBhabhaCrossSection::~G4PolarizedBhabhaCrossSection()
{
}
void G4PolarizedBhabhaCrossSection::Initialize(
			  G4double e,
			  G4double gamma,
			  G4double phi,
		    const G4StokesVector & pol0,
                    const G4StokesVector & pol1,
			  int flag)
{
  G4double re2 = classic_electr_radius * classic_electr_radius;
  G4double gamma2 = gamma*gamma;
  G4double gamma3 = gamma2*gamma;
  G4double gamma4 = gamma2*gamma2;
  G4double gmo = (gamma - 1.);
  G4double gmo2 = (gamma - 1.)*(gamma - 1.);
  G4double gmo3 = gmo2*(gamma - 1.);
  G4double gpo = (gamma + 1.);
  G4double gpo2 = (gamma + 1.)*(gamma + 1.);
  G4double gpo3 = gpo2*(gamma + 1.);
  G4double gpo12 = std::sqrt(gpo);
  G4double gpo32 = gpo*gpo12;
  G4double gpo52 = gpo2*gpo12;

  G4double pref = re2/(gamma - 1.0);
  G4double sqrttwo=std::sqrt(2.);
  G4double d = std::sqrt(1./e - 1.);
  G4double e2 = e*e;
  G4double e3 = e2*e;
  G4double e4 = e2*e2;

  G4double coeff = 0.;

  G4bool polarized=(!pol0.IsZero())||(!pol1.IsZero());

  if (flag==0) polarized=false;
  // Unpolarised part of XS
  phi0 = 0.;
  phi0+= e2*gmo3/gpo3;
  phi0+= -2.*e*gamma*gmo2/gpo3;
  phi0+= (3.*gamma2 + 6.*gamma + 4.)*gmo/gpo3;
  phi0+= -(2.*gamma2 + 4.*gamma + 1.)/(e*gpo2);
  phi0+= gamma2/(e2*(gamma2 - 1.));
  phi0*=0.25;
  // Initial state polarisarion dependence
  if (polarized) {
    G4double usephi=1.;
    if (flag<=1) usephi=0.;
    //    G4cout<<"Polarized differential Bhabha cross section"<<G4endl;
    //    G4cout<<"Initial state polarisation contributions"<<G4endl;
    //    G4cout<<"Diagonal Matrix Elements"<<G4endl;
    G4double xx = 0.;
	xx += -usephi*std::cos(2.*phi)*(e - 1.)*(2.*e2*gmo2 + gamma - 2.*e*(gamma2 - 1.)+ 1.);
	xx += -(4.*gmo2*e3 - 8.*gmo*gamma*e2 + (gamma2 + 3.)*e + 3.*gamma2 + 4.*gamma + 1.)/gpo;
	xx /= 8.*e*gpo2;
    G4double yy = 0.;
	yy += usephi*std::cos(2.*phi)*(e - 1.)*(2.*e2*gmo2 + gamma - 2.*e*(gamma2 - 1.)+ 1.);
	yy += -(4.*gmo2*e3 - 8.*gmo*gamma*e2 + (gamma2 + 3.)*e + 3.*gamma2 + 4.*gamma + 1.)/gpo;
	yy /= 8.*e*gpo2;
    G4double zz = 0.;
	zz += gmo2*(gamma + 3.)*e3 - 2.*gamma*(gamma2 + 2.*gamma - 3.)*e2;
	zz += (3.*gamma3 + 4.*gamma2 - 2.*gamma - 1.)*e - gamma*(2.*gamma2 + 3.*gamma + 1.);
	zz /= 4.*e*gpo3;
    phi0 += xx*pol0.x()*pol1.x() + yy*pol0.y()*pol1.y() + zz*pol0.z()*pol1.z();

    if (usephi==1.) {
    //    G4cout<<"Non-diagonal Matrix Elements"<<G4endl;
      G4double xy, yx;
	xy = (1. - e)*(2.*e2*gmo2 + gamma - 2.*e*(gamma2 - 1.) + 1.)*std::sin(2.*phi);
	xy /= 8.*e*gpo2;
	yx = xy;

      G4double xz, zx, zy, yz;
		coeff = d*(2.*e2*gmo2 + gamma*gpo + e*(-3.*gamma2 + 2.*gamma + 1.));
		coeff /= 2.*sqrttwo*gpo52;
	xz = zx = coeff*std::cos(phi);
	yz = zy = coeff*std::sin(phi);
      phi0+=yx*pol0.y()*pol1.x() + xy*pol0.x()*pol1.y();
      phi0+=zx*pol0.z()*pol1.x() + xz*pol0.x()*pol1.z();
      phi0+=zy*pol0.z()*pol1.y() + yz*pol0.y()*pol1.z();
    }
  }
  // Final state polarisarion dependence
  phi2=G4ThreeVector();
  phi3=G4ThreeVector();

  if (flag>=1) {
    G4double w = (e*gmo - gpo);
    G4double q = (e*gmo + 2.);
    //
    // Final Positron Ppl
    //
	// initial positron Kpl
    if (!pol0.IsZero()) {
      G4double xxPplKpl = gmo2*(4.*gmo*e3 + (3. - 7.*gamma)*e2 + 2.*gamma*e + gamma + 1.)* e*std::cos(2.*phi);
	       xxPplKpl += 4.*gmo3*e4 + gmo2*(2.*gamma2 - 3.*gamma - 3.)*e3;
	       xxPplKpl += -2.*gmo2*(3.*gamma2 + 7.*gamma + 4.)*e2;
	       xxPplKpl += (6.*gamma4 + 9.*gamma3 - 5.*gamma2 - 9.*gamma - 1.)*e;
	       xxPplKpl += -2.*gamma2*gpo2;
	       xxPplKpl /= 8.*e2*gmo*gpo2*w;
      G4double yyPplKpl = -gmo2*(4.*gmo*e3 + (3. - 7.*gamma)*e2 + 2.*gamma*e + gamma + 1.)* e*std::cos(2.*phi);
	       yyPplKpl += 4.*gmo3*e4 + gmo2*(2.*gamma2 - 3.*gamma - 3.)*e3;
	       yyPplKpl += -2.*gmo2*(3.*gamma2 + 7.*gamma + 4.)*e2;
	       yyPplKpl += (6.*gamma4 + 9.*gamma3 - 5.*gamma2 - 9.*gamma - 1.)*e;
	       yyPplKpl += -2.*gamma2*gpo2;
	       yyPplKpl /= 8.*e2*gmo*gpo2*w;
      G4double zzPplKpl = -2.*gmo3*gpo*e4 + gmo2*(5.*gamma2 + gamma - 2.)*e3;
	       zzPplKpl += -gmo2*(5.*gamma2 + 8.*gamma + 4.)*e2;
	       zzPplKpl += gamma*(3.*gamma3 + 5.*gamma2 - 3.*gamma - 5.)*e;
	       zzPplKpl += -gamma2*gpo2;
	       zzPplKpl /= 4.*e2*gmo*gpo2*w;
      G4double xyPplKpl = (e - 1.)*gmo*(4.*gmo*e2 - (3.*gamma + 1.)*e - gamma - 1.)*std::sin(2.*phi);
	       xyPplKpl /= 8.*e*gpo2*w;
      G4double yxPplKpl = xyPplKpl;
	       coeff = d*(-2.*gmo2*e3 + (gamma2 - 1.)*e2 - gamma*e + e + gamma*gpo);
	       coeff /= 2.*sqrttwo*e*gpo32*w;
      G4double xzPplKpl = coeff*std::cos(phi);
      G4double yzPplKpl = coeff*std::sin(phi);
	       coeff = -d*(2.*gmo2*e3 + (-3.*gamma2 + 4.*gamma - 1.)*e2 - gamma*e + e + gamma*gpo);
	       coeff /= 2.*sqrttwo*e*gpo32*w;
      G4double zxPplKpl = coeff*std::cos(phi);
      G4double zyPplKpl = coeff*std::sin(phi);
      phi2[0] += xxPplKpl*pol0.x() + xyPplKpl*pol0.y() + xzPplKpl*pol0.z();
      phi2[1] += yxPplKpl*pol0.x() + yyPplKpl*pol0.y() + yzPplKpl*pol0.z();
      phi2[2] += zxPplKpl*pol0.x() + zyPplKpl*pol0.y() + zzPplKpl*pol0.z();
    }
	// initial electron Kmn
    if (!pol1.IsZero()) {
      G4double xxPplKmn = 4.*gmo2*e3 - (2.*gamma3 + 7.*gamma2 - 12.*gamma + 3.)*e2;
	       xxPplKmn += 2.*gamma*(gamma2 + 2.*gamma - 1.)*e - gpo2;
	       xxPplKmn *= std::cos(2.*phi);
	       xxPplKmn += 4.*gmo2*e3 + (-5.*gamma2 + 2.*gamma + 3.)*e2;
	       xxPplKmn += (-6.*gamma2 + 2.*gamma + 8.)*e + 3.*gamma2 + 4.*gamma + 1.;
	       xxPplKmn /= 8.*e*gpo2*w;
      G4double yyPplKmn = 4.*gmo2*e3 - (2.*gamma3 + 7.*gamma2 - 12.*gamma + 3.)*e2;
	       yyPplKmn += 2.*gamma*(gamma2 + 2.*gamma - 1.)*e - gpo2;
	       yyPplKmn *= -std::cos(2.*phi);
	       yyPplKmn += 4.*gmo2*e3 + (-5.*gamma2 + 2.*gamma + 3.)*e2;
	       yyPplKmn += (-6.*gamma2 + 2.*gamma + 8.)*e + 3.*gamma2 + 4.*gamma + 1.;
	       yyPplKmn /= 8.*e*gpo2*w;
      G4double zzPplKmn = -2.*gmo2*gpo*e3 + (5.*gamma3 - 7.*gamma + 2.)*e2;
	       zzPplKmn += (-5.*gamma3 - 7.*gamma2 + 4.*gpo)*e + gamma*(2.*gamma2 + 3.*gamma + 1.);
	       zzPplKmn /= 4.*e*gpo2*w;
      G4double xyPplKmn = 4.*gmo2*e3 - (2.*gamma3 + 7.*gamma2 - 12.*gamma + 3.)*e2;
	       xyPplKmn += 2.*gamma*(gamma2 + 2.*gamma - 1.)*e - gpo2;
	       xyPplKmn *= std::sin(2.*phi);
	       xyPplKmn /= 8.*e*gpo2*w;
      G4double yxPplKmn = xyPplKmn;
	       coeff = -d*(2.*e2*gmo2 + e*(-3.*gamma2 + 2.*gamma + 1.)+ 3.*gamma2 - 1.);
	       coeff /= 2.*sqrttwo*gpo32*w;
      G4double xzPplKmn = coeff*std::cos(phi);
      G4double yzPplKmn = coeff*std::sin(phi);
	       coeff = -d*(2.*e2*gmo2 + e*(-5.*gamma2 + 6.*gamma - 1.) + gamma2 + 1.);
	       coeff /= 2.*sqrttwo*gpo32*w;
      G4double zxPplKmn = coeff*std::cos(phi);
      G4double zyPplKmn = coeff*std::sin(phi);
      phi2[0] += xxPplKmn*pol1.x() + xyPplKmn*pol1.y() + xzPplKmn*pol1.z();
      phi2[1] += yxPplKmn*pol1.x() + yyPplKmn*pol1.y() + yzPplKmn*pol1.z();
      phi2[2] += zxPplKmn*pol1.x() + zyPplKmn*pol1.y() + zzPplKmn*pol1.z();
    }
//
// Final Electron Pmn
//
	// initial positron Kpl
    if (!pol0.IsZero()) {
      G4double xxPmnKpl = -(e - 1.)*(4.*e2*gmo2 + e*(-5.*gamma2 + 2.*gamma + 3.) + 2.*gpo)*std::cos(2.*phi);
	       xxPmnKpl += -4.*gmo2*e3 + (2.*gamma3 + 11.*gamma2 - 20.*gamma + 7.)*e2;
	       xxPmnKpl += (-5.*gamma2 + 10.*gamma - 9.)*e - 6.*gamma - 2.;
	       xxPmnKpl /= 8.*e*gpo2*q;
      G4double yyPmnKpl = (e - 1.)*(4.*e2*gmo2 + e*(-5.*gamma2 + 2.*gamma + 3.) + 2.*gpo)*std::cos(2.*phi);
	       yyPmnKpl += -4.*gmo2*e3 + (2.*gamma3 + 11.*gamma2 - 20.*gamma + 7.)*e2;
	       yyPmnKpl += (-5.*gamma2 + 10.*gamma - 9.)*e - 6.*gamma - 2.;
	       yyPmnKpl /= 8.*e*gpo2*q;
      G4double zzPmnKpl = 2.*gmo2*gpo*e3 + gamma*(-3.*gamma2 - 2.*gamma + 5.)*e2;
	       zzPmnKpl += (2.*gamma3 + 9.*gamma2 - 3.*gamma - 4.)*e - 2.*gamma*(2.*gamma + 1.);
	       zzPmnKpl /= 4.*e*gpo2*q;
      G4double xyPmnKpl = -(e - 1.)*(4.*e2*gmo2 + (-5.*gamma2 + 2.*gamma + 3.)*e + 2.*gpo)*std::sin(2.*phi);
	       xyPmnKpl /= 8.*e*gpo2*q;
      G4double yxPmnKpl = xyPmnKpl;
	       coeff = d*(2.*e2*gmo2 + (-5.*gamma2 + 6.*gamma - 1.)*e + 4.*gamma2 - 2.);
	       coeff /= 2.*sqrttwo*gpo32*q;
      G4double xzPmnKpl = coeff*std::cos(phi);
      G4double yzPmnKpl = coeff*std::sin(phi);
	       coeff = d*(2.*e2*gmo2 + (-3.*gamma2 + 2.*gamma + 1.)*e + 2.);
	       coeff /= 2.*sqrttwo*gpo32*q;
      G4double zxPmnKpl = coeff*std::cos(phi);
      G4double zyPmnKpl = coeff*std::sin(phi);
      phi3[0] += xxPmnKpl*pol0.x() + xyPmnKpl*pol0.y() + xzPmnKpl*pol0.z();
      phi3[1] += yxPmnKpl*pol0.x() + yyPmnKpl*pol0.y() + yzPmnKpl*pol0.z();
      phi3[2] += zxPmnKpl*pol0.x() + zyPmnKpl*pol0.y() + zzPmnKpl*pol0.z();
    }
	// initial electron Kmn
    if (!pol1.IsZero()) {
      G4double xxPmnKmn = -(e - 1.)*gmo2*(4.*gmo*e2 + (2.*gamma2 - 3.*gpo)*e - 2.*gpo2)*e*std::cos(2.*phi);
	       xxPmnKmn += -4.*gmo3*e4 + gmo2*(11.*gamma - 7.)*e3 - 3.*(gamma - 3.)*gmo2*e2;
	       xxPmnKmn += (-8.*gamma3 - 6.*gamma2 +12.*gamma + 2.)*e +4.*gamma2*gpo;
	       xxPmnKmn /= 8.*e2*gmo*gpo2*q;
      G4double yyPmnKmn = (e - 1.)*gmo2*(4.*gmo*e2 + (2.*gamma2 - 3.*gpo)*e - 2.*gpo2)*e*std::cos(2.*phi);
	       yyPmnKmn += -4.*gmo3*e4 + gmo2*(11.*gamma - 7.)*e3 - 3.*(gamma - 3.)*gmo2*e2;
	       yyPmnKmn += (-8.*gamma3 - 6.*gamma2 +12.*gamma + 2.)*e +4.*gamma2*gpo;
	       yyPmnKmn /= 8.*e2*gmo*gpo2*q;
      G4double zzPmnKmn = 2.*gmo3*gpo*e4 - gmo2*gamma*(3.*gamma + 1.)*e3;
	       zzPmnKmn += gmo2*(2.*gamma2 + 3.*gamma + 4.)*e2;
	       zzPmnKmn += -gamma*(gamma3 + 4.*gamma2 + gamma - 6.)*e +2.*gamma2*gpo;
	       zzPmnKmn /= 4.*e2*gmo*gpo2*q;
      G4double xyPmnKmn = -(e - 1.)*gmo*std::sin(2.*phi);
	       xyPmnKmn *= 4.*gmo*e2 + (2.*gamma2 - 3.*gpo)*e - 2.*gpo2;
	       xyPmnKmn /= 8.*e*gpo2*q;
      G4double yxPmnKmn = xyPmnKmn;
	       coeff = d*(2.*gmo2*e3 + (-3.*gamma2 + 4.*gamma - 1.)*e2 
			  + (3.*gamma2 - gamma - 2.)*e - 2.*gamma*gpo);
	       coeff /= 2.*sqrttwo*e*gpo32*q;
      G4double xzPmnKmn = coeff*std::cos(phi);
      G4double yzPmnKmn = coeff*std::sin(phi);
	       coeff = d*(2.*gmo2*e3 - (gamma2 - 1.)*e2 + (-3.*gamma2 + gamma + 2.)*e + 2.*gamma*gpo);
	       coeff /= 2.*sqrttwo*e*gpo32*q;
      G4double zxPmnKmn = coeff*std::cos(phi);
      G4double zyPmnKmn = coeff*std::sin(phi);
      phi3[0] += xxPmnKmn*pol1.x() + xyPmnKmn*pol1.y() + xzPmnKmn*pol1.z();
      phi3[1] += yxPmnKmn*pol1.x() + yyPmnKmn*pol1.y() + yzPmnKmn*pol1.z();
      phi3[2] += zxPmnKmn*pol1.x() + zyPmnKmn*pol1.y() + zzPmnKmn*pol1.z();
    }
  }
phi0 *= pref;
phi2 *= pref;
phi3 *= pref;

}

G4double G4PolarizedBhabhaCrossSection::XSection(const G4StokesVector & pol2,const G4StokesVector & pol3)
{
  G4double xs=0.;
  xs+=phi0;

  G4bool polarized=(!pol2.IsZero())||(!pol3.IsZero());
  if (polarized) {
    xs+=phi2*pol2 + phi3*pol3;
  }
  return xs;
}

G4double G4PolarizedBhabhaCrossSection::TotalXSection(
  const G4double xmin, const G4double xmax, const G4double gamma,
  const G4StokesVector & pol0,const G4StokesVector & pol1)
{
  G4double xs=0.;

  G4double x=xmin;

  if (xmax != 1.) G4cout<<" warning xmax expected to be 1 but is "<<xmax<< G4endl;

  //  re -> electron radius^2;
  G4double re2 = classic_electr_radius * classic_electr_radius;
  G4double gamma2=gamma*gamma;
  G4double gmo2 = (gamma - 1.)*(gamma - 1.);
  G4double gpo2 = (gamma + 1.)*(gamma + 1.);
  G4double gpo3 = gpo2*(gamma + 1.);
  G4double logMEM = std::log(x);
  G4double pref = twopi*re2/(gamma - 1.0);
  // unpolarise XS
  G4double sigma0 = 0.;
	sigma0 += -gmo2*(gamma - 1.)*x*x*x/3. + gmo2*gamma*x*x;
	sigma0 += -(gamma - 1.)*(3.*gamma*(gamma + 2.) +4.)*x;
	sigma0 += (gamma*(gamma*(gamma*(4.*gamma - 1.) - 21.) - 7.)+13.)/(3.*(gamma - 1.));
	sigma0 /= gpo3;
	sigma0 += logMEM*(2. - 1./gpo2);
	sigma0 += gamma2/((gamma2 - 1.)*x);
  //    longitudinal part
  G4double sigma2=0.;
  	sigma2 += logMEM*gamma*(gamma + 1.)*(2.*gamma + 1.);
	sigma2 += gamma*(7.*gamma*(gamma + 1.) - 2.)/3.;
	sigma2 += -(3.*gamma + 1.)*(gamma2 + gamma - 1.)*x;
	sigma2 += (gamma - 1.)*gamma*(gamma + 3.)*x*x;
	sigma2 += -gmo2*(gamma + 3.)*x*x*x/3.;
	sigma2 /= gpo3;
  //    transverse part
  G4double sigma3=0.;
  	sigma3 += 0.5*(gamma + 1.)*(3.*gamma + 1.)*logMEM;
  	sigma3 += (gamma*(5.*gamma - 4.) - 13.)/6.;
  	sigma3 += 0.5*(gamma2 + 3.)*x;
  	sigma3 += - 2.*(gamma - 1.)*gamma*x*x;  // *AS* changed sign
  	sigma3 += 2.*gmo2*x*x*x/3.;
	sigma3 /= gpo3;
  // total cross section
  xs+=pref*(sigma0 + sigma2*pol0.z()*pol1.z() + sigma3*(pol0.x()*pol1.x()+pol0.y()*pol1.y()));

  return xs;
}


G4StokesVector G4PolarizedBhabhaCrossSection::GetPol2()
{
  // Note, mean polarization can not contain correlation
  // effects.
  return  1./phi0 * phi2;
}
G4StokesVector G4PolarizedBhabhaCrossSection::GetPol3()
{
  // Note, mean polarization can not contain correlation
  // effects.
  return  1./phi0 * phi3;
}
