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
// $Id: G4VPolarizedCrossSection.cc 91742 2015-08-04 11:48:51Z gcosmo $
// File name:     G4VPolarizedCrossSection
//
// Author:        Andreas Schaelicke
//
// Creation date: 15.05.2005
//
// Modifications:
//
// Class Description:
//   (pure virtual) interface class
//
//   provides readable but efficient routines to determine 
//   polarization for the final state of a given process
//   empoying the differential cross section
//

#include "G4VPolarizedCrossSection.hh"
#include "Randomize.hh"

G4VPolarizedCrossSection::G4VPolarizedCrossSection() :
  fXmin(0), fXmax(1.), fYmin(1.), theA(1), theZ(1), fCoul(0.)
{
}

G4VPolarizedCrossSection::~G4VPolarizedCrossSection()
{
}

void G4VPolarizedCrossSection::Initialize(G4double, G4double, G4double,
					  const G4StokesVector &,
					  const G4StokesVector &,
					  G4int ) 
{
} 
 
G4StokesVector G4VPolarizedCrossSection::GetPol2()
{
  // neglects correlation effects! 

  G4double invXsecTotal=1./XSection(G4StokesVector::ZERO,G4StokesVector::ZERO);
  G4double xsPol1=XSection(G4StokesVector::P1,G4StokesVector::ZERO);
  G4double xsPol2=XSection(G4StokesVector::P2,G4StokesVector::ZERO);
  G4double xsPol3=XSection(G4StokesVector::P3,G4StokesVector::ZERO);
  return G4ThreeVector(invXsecTotal*xsPol1,invXsecTotal*xsPol2,invXsecTotal*xsPol3);
}

G4StokesVector G4VPolarizedCrossSection::GetPol3()
{
  // neglects correlation effects! 

  G4double invXsecTotal=1./XSection(G4StokesVector::ZERO,G4StokesVector::ZERO);
  G4double xsPol1=XSection(G4StokesVector::ZERO,G4StokesVector::P1); 
  G4double xsPol2=XSection(G4StokesVector::ZERO,G4StokesVector::P2); 
  G4double xsPol3=XSection(G4StokesVector::ZERO,G4StokesVector::P3); 
  return G4ThreeVector(invXsecTotal*xsPol1,invXsecTotal*xsPol2,invXsecTotal*xsPol3); 
}
// minimal energy fraction in TotalXSection
G4double G4VPolarizedCrossSection::GetXmin(G4double /*y*/)
{
  return fXmin;
}

// maximal energy fraction in TotalXSection
G4double G4VPolarizedCrossSection::GetXmax(G4double /*y*/)
{
  return fXmax;
}


/*
void G4VPolarizedCrossSection::DicePolarization() 
{
  // can respect correlation effects, but is limited to 
  // one quantization axis!
  G4double sigma[4];
  sigma[0]=XSection(G4StokesVector::P3,G4StokesVector::P3);
  sigma[1]=XSection(G4StokesVector::P3,G4StokesVector::M3);
  sigma[2]=XSection(G4StokesVector::M3,G4StokesVector::P3);
  sigma[3]=XSection(G4StokesVector::M3,G4StokesVector::M3);

  G4double sigma_max = 4. * XSection(G4StokesVector::ZERO,G4StokesVector::ZERO);

  for (G4int i=0;i<4;++i) {
    G4cout<<"sigma="<<sigma[i]<<" vs."<<(.25*sigma_max)<<G4endl;
    if (sigma[i]<0 || sigma[i]>sigma_max) {
      G4cout<<"ERROR G4VPolarizedCrossSection::DicePolarization(["<<i<<"]):  "
	    <<sigma[i]<<" vs."<<sigma_max<<G4endl;
    }
    if (i>0) sigma[i]+=sigma[i-1];
  }

  G4int k = 0;
  G4double disc = sigma[3]*G4UniformRand();
  // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  while (sigma[k]<disc && k<4) {
    ++k;
  }

  if ((k&2)==0) pol2=G4StokesVector::P3;
  else          pol2=G4StokesVector::M3;
  if ((k&1)==0) pol3=G4StokesVector::P3;
  else          pol3=G4StokesVector::M3;

}
*/

/*
G4StokesVector G4VPolarizedCrossSection::DicedPol2()
{
  return pol2;
}

G4StokesVector G4VPolarizedCrossSection::DicedPol3()
{
  return pol3;
}
*/

G4double G4VPolarizedCrossSection::TotalXSection(G4double, G4double, G4double,
						 const G4StokesVector &,const G4StokesVector &)
{
  G4cout << "WARNING virtual function G4VPolarizedCrossSection::TotalXSection() called" << G4endl;
  return 0.;
}
