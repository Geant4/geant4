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


#include "G4NeutrinoElectronTotXsc.hh"
#include "G4NeutrinoElectronCcXsc.hh"
#include "G4NeutrinoElectronNcXsc.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"
using namespace std;
using namespace CLHEP;

G4NeutrinoElectronTotXsc::G4NeutrinoElectronTotXsc()
 : G4VCrossSectionDataSet("NuElectronTotXsc")
{
  fCcXsc = new G4NeutrinoElectronCcXsc();
  fNcXsc = new G4NeutrinoElectronNcXsc();

  fCutEnergy = 0.; // default value
  fBiasingFactor = 1.; // default as physics
  fCcRatio = 0.5;
}

G4NeutrinoElectronTotXsc::~G4NeutrinoElectronTotXsc() 
{}

//////////////////////////////////////////////////////

G4bool 
G4NeutrinoElectronTotXsc::IsElementApplicable( const G4DynamicParticle* aPart, G4int i, const G4Material* mat)
{
  G4bool result  = false;

  G4bool apCc = fCcXsc->IsElementApplicable( aPart,i, mat);
  G4bool apNc = fNcXsc->IsElementApplicable( aPart, i, mat);

  if( apCc || apNc) result = true;

  return result;
}

////////////////////////////////////////////////////

G4double G4NeutrinoElectronTotXsc::
GetElementCrossSection(const G4DynamicParticle* aPart, G4int ZZ,  
		       const G4Material* mat) 
{
  G4double result = 0.;

  G4double ccXsc = fCcXsc->GetElementCrossSection( aPart, ZZ, mat); 
  G4double ncXsc = fNcXsc->GetElementCrossSection( aPart, ZZ, mat); 

  result = ccXsc + ncXsc;

  if (result > 0.) fCcRatio = ccXsc/result;
  else             fCcRatio = 0.;

  return result;
}

////////////////////////////////////////////////////

void G4NeutrinoElectronTotXsc::SetBiasingFactor(G4double bf)
{
    fBiasingFactor = bf;
    fCcXsc->SetBiasingFactor(bf);
    fNcXsc->SetBiasingFactor(bf);
}

////////////////////////////////////////////////////
//
// For separate testing Cc and Nc currents

void G4NeutrinoElectronTotXsc::SetBiasingFactors(G4double bfCc, G4double bfNc)
{
    fCcXsc->SetBiasingFactor(bfCc);
    fNcXsc->SetBiasingFactor(bfNc);
}
