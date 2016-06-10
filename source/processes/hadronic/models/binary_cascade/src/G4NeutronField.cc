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
//      GEANT 4 class implementation file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4NeutronField.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#include "G4NeutronField.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VNuclearDensity.hh"
#include "G4FermiMomentum.hh"

G4NeutronField::G4NeutronField(G4V3DNucleus * aNucleus) : 
    G4VNuclearField(aNucleus), theDensity(theNucleus->GetNuclearDensity())
{ 
  theA = theNucleus->GetMassNumber();
  theZ = theNucleus->GetCharge();
  theFermi.Init(theA, theZ);
  theR = 2.*theNucleus->GetOuterRadius();
  for (G4double aR=0.;aR<theR; aR+=0.3*fermi)
  {
    G4ThreeVector aPosition(0,0,aR);
    G4double density = GetDensity(aPosition);
    G4double fermiMom = GetFermiMomentum(density);
    theFermiMomBuffer.push_back(fermiMom);
  }
  {
  G4ThreeVector aPosition(0,0,theR);
  G4double density = GetDensity(aPosition);
  G4double fermiMom = GetFermiMomentum(density);
  theFermiMomBuffer.push_back(fermiMom);  
  }
  {
  G4ThreeVector aPosition(0,0,theR+0.001*fermi);
  theFermiMomBuffer.push_back(0);  
  }
  {
  G4ThreeVector aPosition(0,0,1.*m);
  theFermiMomBuffer.push_back(0);  
  }
}

G4NeutronField::~G4NeutronField()
{ }

G4double G4NeutronField::GetField(const G4ThreeVector & aPosition)
{
  G4double x = aPosition.mag();
  unsigned int index = static_cast<unsigned int>(x/(0.3*fermi));
  if( (index+2) > theFermiMomBuffer.size()) return theFermiMomBuffer.back();
  G4double y1 = theFermiMomBuffer[index];
  G4double y2 = theFermiMomBuffer[index+1];
  G4double x1 = (0.3*fermi)*index;
  G4double x2 = (0.3*fermi)*(index+1);
  G4double fermiMom = y1 + (x-x1)*(y2-y1)/(x2-x1);
  return -1*(fermiMom*fermiMom)/(2*neutron_mass_c2);
}

G4double G4NeutronField::GetBarrier()
{
/*
 *   G4double A = theNucleus->GetMassNumber();
 *   G4double Z = theNucleus->GetCharge();
 * 
 *   return G4NucleiPropertiesTable::GetBindingEnergy(Z, A)/A;
 */
	return 0.;
}






