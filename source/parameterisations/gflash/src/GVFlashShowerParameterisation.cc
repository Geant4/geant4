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
//
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ------- GVFlashShowerParameterisation -------
//
// Authors: Joanna Weng - 11.2005
// ------------------------------------------------------------

#include <cmath>

#include "GVFlashShowerParameterisation.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Material.hh"
#include "Gamma.hh" // @@@@
#include "G4MaterialTable.hh"

GVFlashShowerParameterisation::GVFlashShowerParameterisation()
  : thePar(0), density(0.), A(0.), Z(0.), X0(0.), Ec(0.), Rm(0.), NSpot(0.)
{
  fGamma = new MyGamma;
}

GVFlashShowerParameterisation::~GVFlashShowerParameterisation()
{
  delete fGamma;
}

G4double GVFlashShowerParameterisation::GetEffZ(const G4Material * mat  )
{
  // Returns Z or effective Z=sum(pi*Zi) (if compound/mixture)
  // of given material
  //
  G4double z = 0.;
  G4int nofElements = (G4int)mat->GetNumberOfElements();
  if (nofElements > 1) 
  {
    for (G4int i=0; i<nofElements; ++i) {
      G4double zOfElement = mat->GetElement(i)->GetZ();
      G4double massFraction = mat->GetFractionVector()[i];
      // cout << mat->GetElement(i)->GetName()
      //      <<" Z= "<<zOfElement << " , Fraction= "<<massFraction <<endl;
      z += zOfElement*massFraction;
    }
  }
  else { 
    z = mat->GetZ(); 
  }  
  return z;
}

G4double GVFlashShowerParameterisation::GetEffA  (const G4Material * mat  )
{
  // Returns A or effective A=sum(pi*Ai) (if compound/mixture)
  // of given material
  //
  G4double a = 0.;
  G4int nofElements = (G4int)mat->GetNumberOfElements();
  if (nofElements > 1) {
    for (G4int i=0; i<nofElements; ++i) {
      G4double aOfElement = mat->GetElement(i)->GetA()/(g/mole);
      G4double massFraction = mat->GetFractionVector()[i];     
      a += aOfElement*massFraction;
    }
  }
  else { 
    a = mat->GetA()/(g/mole);
  }
  return a;
}

void GVFlashShowerParameterisation::PrintMaterial(const G4Material * mat)
{
  G4cout<<"/********************************************/ " << G4endl;
  G4cout<<"  - GVFlashShowerParameterisation::Material -  " << G4endl;
  G4cout<<"        Material : " << mat->GetName()  << G4endl;
  G4cout<<"   Z  = " << Z      << G4endl;
  G4cout<<"   A  = " << A      << G4endl;
  G4cout<<"   X0 = " << X0/cm  << " cm" << G4endl;
  G4cout<<"   Rm = " << Rm/cm  << " cm" << G4endl;
  G4cout<<"   Ec = " << Ec/MeV << " MeV"<< G4endl;
  G4cout<<"/********************************************/ " << G4endl; 
}

G4double GVFlashShowerParameterisation::GeneratePhi()
{
  G4double Phi = twopi*G4UniformRand() ;
  return Phi;
}

G4double GVFlashShowerParameterisation::gam(G4double x, G4double a) const 
{
  return fGamma->Gamma(a, x); 
}
