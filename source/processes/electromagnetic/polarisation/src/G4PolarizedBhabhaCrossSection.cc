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
// $Id: G4PolarizedBhabhaCrossSection.cc 68046 2013-03-13 14:31:38Z gcosmo $
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
//   07-11-06 modify reference system for polarisation vectors 
//            (A. Schaelicke & P.Starovoitov)
//
// Class Description:
//   * calculates the differential cross section
//     incomming positron Kpl(along positive z direction) scatters at 
//     an electron Kmn at rest
//   * phi denotes the angle between the scattering plane (defined by the
//     outgoing electron) and X-axis
//   * all stokes vectors refer to spins in the Global System (X,Y,Z)
//

#include "G4PolarizedBhabhaCrossSection.hh"
#include "G4PhysicalConstants.hh"

G4PolarizedBhabhaCrossSection::G4PolarizedBhabhaCrossSection() : phi0(1.)
{
}
G4PolarizedBhabhaCrossSection::~G4PolarizedBhabhaCrossSection()
{
}
void G4PolarizedBhabhaCrossSection::Initialize(
			  G4double e,
			  G4double gamma,
			  G4double /*phi*/,
		    const G4StokesVector & pol0,
                    const G4StokesVector & pol1,
			  G4int flag)
{
  SetXmax(1.);

  G4double re2 = classic_electr_radius * classic_electr_radius;
  G4double gamma2 = gamma*gamma;
  G4double gamma3 = gamma2*gamma;
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

  // *** new ***
  G4double gmo12 = std::sqrt(gmo);
  G4double gmo32 = gmo*gmo12;
  G4double egmp32 = std::pow(e*(2 + e*gmo)*gpo,(3./2.));
  G4double e32 = e*std::sqrt(e);

  G4bool polarized=(!pol0.IsZero())||(!pol1.IsZero());

  if (flag==0) polarized=false;
  // Unpolarised part of XS
  // *AS* UnpME . OK 
  phi0 = 0.;
  phi0+= e2*gmo3/gpo3;
  phi0+= -2.*e*gamma*gmo2/gpo3;
  phi0+= (3.*gamma2 + 6.*gamma + 4.)*gmo/gpo3;
  phi0+= -(2.*gamma2 + 4.*gamma + 1.)/(e*gpo2);
  phi0+= gamma2/(e2*(gamma2 - 1.));
  phi0*=0.25;
  // Initial state polarisarion dependence
  if (polarized) {
    //    G4cout<<"Polarized differential Bhabha cross section"<<G4endl;
    //    G4cout<<"Initial state polarisation contributions"<<G4endl;
    //    G4cout<<"Diagonal Matrix Elements"<<G4endl;
    // *** new ***
    G4double xx = -((e*gmo - gamma)*(-1 - gamma + e*(e*gmo - gamma)*(3 + gamma)))/(4*e*gpo3);
    G4double yy = (e3*gmo3 - 2*e2*gmo2*gamma - gpo*(1 + 2*gamma) + e*(-2 + gamma2 + gamma3))/(4*e*gpo3);
    G4double zz = ((e*gmo - gamma)*(e2*gmo*(3 + gamma) - e*gamma*(3 + gamma) + gpo*(1 + 2*gamma)))/(4*e*gpo3);
    // ***

    phi0 += xx*pol0.x()*pol1.x() + yy*pol0.y()*pol1.y() + zz*pol0.z()*pol1.z();

    {
      G4double xy = 0;
      G4double xz = (d*(e*gmo - gamma)*(-1 + 2*e*gmo - gamma))/(2*sqrttwo*gpo52);
      G4double yx = 0;
      G4double yz = 0;
      G4double zx = xz;
      G4double zy = 0;
    //    G4cout<<"Non-diagonal Matrix Elements"<<G4endl;
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
    // Final Positron Ppl
    //
	// initial positron Kpl
    if (!pol0.IsZero()) {

      G4double xxPplKpl = -((-1 + e)*(e*gmo - gamma)*(-(gamma*gpo) +  e*(-2 + gamma + gamma2)))/
	(4*e2*gpo*std::sqrt(gmo*gpo*(-1 + e + gamma - e*gamma)*  (1 + e + gamma - e*gamma)));
      G4double xyPplKpl = 0;
      G4double xzPplKpl = ((e*gmo - gamma)*(-1 - gamma + e*gmo*(1 + 2*gamma)))/
	(2*sqrttwo*e32*gmo*gpo2*std::sqrt(1 + e + gamma - e*gamma));
      G4double yxPplKpl = 0;
      G4double yyPplKpl = (gamma2*gpo + e2*gmo2*(3 + gamma) - 
		  e*gmo*(1 + 2*gamma*(2 + gamma)))/(4*e2*gmo*gpo2);
      G4double yzPplKpl = 0;
      G4double zxPplKpl = ((e*gmo - gamma)*(1 + e*(-1 + 2*e*gmo - 2*gamma)*gmo + gamma))/
	(2*sqrttwo*e*gmo*gpo2*std::sqrt(e*(1 + e + gamma - e*gamma)));
      G4double zyPplKpl = 0;
      G4double zzPplKpl = -((e*gmo - gamma)*std::sqrt((1 - e)/(e - e*gamma2 + gpo2))*
		   (2*e2*gmo2 + gamma + gamma2 - e*(-2 + gamma + gamma2)))/
	(4*e2*(-1 + gamma2));

      phi2[0] += xxPplKpl*pol0.x() + xyPplKpl*pol0.y() + xzPplKpl*pol0.z();
      phi2[1] += yxPplKpl*pol0.x() + yyPplKpl*pol0.y() + yzPplKpl*pol0.z();
      phi2[2] += zxPplKpl*pol0.x() + zyPplKpl*pol0.y() + zzPplKpl*pol0.z();
    }
	// initial electron Kmn
    if (!pol1.IsZero()) {
      G4double xxPplKmn = ((-1 + e)*(e*(-2 + gamma)*gmo + gamma))/(4*e*gpo32*std::sqrt(1 + e2*gmo + gamma - 2*e*gamma));
      G4double xyPplKmn = 0;
      G4double xzPplKmn = (-1 + e*gmo + gmo*gamma)/(2*sqrttwo*gpo2* std::sqrt(e*(1 + e + gamma - e*gamma)));
      G4double yxPplKmn = 0;
      G4double yyPplKmn = (-1 - 2*gamma + e*gmo*(3 + gamma))/(4*e*gpo2);
      G4double yzPplKmn = 0;
      G4double zxPplKmn = (1 + 2*e2*gmo2 + gamma + gamma2 + e*(1 + (3 - 4*gamma)*gamma))/
	(2*sqrttwo*gpo2*std::sqrt(e*(1 + e + gamma - e*gamma)));
      G4double zyPplKmn = 0;
      G4double zzPplKmn = -(std::sqrt((1 - e)/(e - e*gamma2 + gpo2))*
		   (2*e2*gmo2 + gamma + 2*gamma2 + e*(2 + gamma - 3*gamma2)))/(4*e*gpo);

      phi2[0] += xxPplKmn*pol1.x() + xyPplKmn*pol1.y() + xzPplKmn*pol1.z();
      phi2[1] += yxPplKmn*pol1.x() + yyPplKmn*pol1.y() + yzPplKmn*pol1.z();
      phi2[2] += zxPplKmn*pol1.x() + zyPplKmn*pol1.y() + zzPplKmn*pol1.z();
    }
//
// Final Electron Pmn
//
	// initial positron Kpl
    if (!pol0.IsZero()) {
      G4double xxPmnKpl = ((-1 + e*gmo)*(2 + gamma))/(4*gpo* std::sqrt(e*(2 + e*gmo)*gpo));
      G4double xyPmnKpl = 0;
      G4double xzPmnKpl = (std::sqrt((-1 + e)/(-2 + e - e*gamma))*
	    (e + gamma + e*gamma -   2*(-1 + e)*gamma2))/(2*sqrttwo*e*gpo2);
      G4double yxPmnKpl = 0;
      G4double yyPmnKpl = (-1 - 2*gamma + e*gmo*(3 + gamma))/(4*e*gpo2);
      G4double yzPmnKpl = 0;
      G4double zxPmnKpl = -((-1 + e)*(1 + 2*e*gmo)*(e*gmo - gamma))/
	(2*sqrttwo*e*std::sqrt(-((-1 + e)*(2 + e*gmo)))*gpo2);
      G4double zyPmnKpl = 0;
      G4double zzPmnKpl = (-2 + 2*e2*gmo2 + gamma*(-1 + 2*gamma) + 
		  e*(-2 + (5 - 3*gamma)*gamma))/(4*std::sqrt(e*(2 + e*gmo))*  gpo32);

      phi3[0] += xxPmnKpl*pol0.x() + xyPmnKpl*pol0.y() + xzPmnKpl*pol0.z();
      phi3[1] += yxPmnKpl*pol0.x() + yyPmnKpl*pol0.y() + yzPmnKpl*pol0.z();
      phi3[2] += zxPmnKpl*pol0.x() + zyPmnKpl*pol0.y() + zzPmnKpl*pol0.z();
    }
	// initial electron Kmn
    if (!pol1.IsZero()) {
      G4double xxPmnKmn = -((2 + e*gmo)*(-1 + e*gmo - gamma)*(e*gmo - gamma)*
		   (-2 + gamma))/(4*gmo*egmp32);
      G4double xyPmnKmn = 0;
      G4double xzPmnKmn = ((e*gmo - gamma)*
			   std::sqrt((-1 + e + gamma - e*gamma)/(2 + e*gmo))*
			   (e + gamma - e*gamma + gamma2))/
	(2*sqrttwo*e2*gmo32*gpo2);
      G4double yxPmnKmn = 0;
      G4double yyPmnKmn = (gamma2*gpo + e2*gmo2*(3 + gamma) - 
		  e*gmo*(1 + 2*gamma*(2 + gamma)))/(4*e2*gmo*gpo2);
      G4double yzPmnKmn = 0;
      G4double zxPmnKmn = -((-1 + e)*(e*gmo - gamma)*(e*gmo + 2*e2*gmo2 - gamma*gpo))/
	(2*sqrttwo*e2*std::sqrt(-((-1 + e)*(2 + e*gmo)))* gmo*gpo2);
      G4double zyPmnKmn = 0;
      G4double zzPmnKmn = ((e*gmo - gamma)*std::sqrt(e/((2 + e*gmo)*gpo))*
		  (-(e*(-2 + gamma)*gmo) + 2*e2*gmo2 + (-2 + gamma)*gpo))/(4*e2*(-1 + gamma2));

      phi3[0] += xxPmnKmn*pol1.x() + xyPmnKmn*pol1.y() + xzPmnKmn*pol1.z();
      phi3[1] += yxPmnKmn*pol1.x() + yyPmnKmn*pol1.y() + yzPmnKmn*pol1.z();
      phi3[2] += zxPmnKmn*pol1.x() + zyPmnKmn*pol1.y() + zzPmnKmn*pol1.z();
    }
  }
  phi0 *= pref;
  phi2 *= pref;
  phi3 *= pref;

}

G4double G4PolarizedBhabhaCrossSection::XSection(const G4StokesVector & pol2,
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

G4double G4PolarizedBhabhaCrossSection::TotalXSection(
  G4double xmin, G4double xmax, G4double gamma,
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
