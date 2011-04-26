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
// $Id: BremLPMTest.cc,v 1.2 2008-08-22 08:16:17 schaelic Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// ------------------------------------------------------------
//   Test Routine for eBremsstrahlung-Models
//      -- including LPM effect
// ------------------------------------------------------------
// run
//   make -f GNUmakefile.root G4TARGET=BremLPMTest
//   BremLPMTest
// this creates a file eBremRel01.root which can be read using e.g.
//   python readBremLPM.py
// ------------------------------------------------------------
//
//  History
//   21.08.08  basic test created (A.Schaelicke)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TTree.h"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"

// ------------------------------------------------------------
// nasty trick to access private or protected functions
// use only in tests and with care !!!
#define  protected public
#define  private public
// ------------------------------------------------------------

#include "G4eBremParametrizedModel.hh"


#include "G4eBremsstrahlungRelModel.hh"
#include "G4eBremsstrahlungModel.hh"


using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void CalcCrossSection()
{
  // initialize models
  G4eBremsstrahlungModel * model1 = new G4eBremsstrahlungModel(); 
  model1->SetLPMFlag(false);
  model1->lowKinEnergy=10.*MeV;
  model1->SetLowEnergyLimit(10.*MeV);

  G4eBremsstrahlungRelModel * modelR = new G4eBremsstrahlungRelModel(); // rel. version
  modelR->SetLPMFlag(false);


  G4eBremParametrizedModel * model2 = new G4eBremParametrizedModel();
  model2->SetLPMFlag(false);
  

  // define energy range and cut
  const G4int nmax=80;
  G4double  emin=1.e2*MeV;
  G4double  emax=1.e6*MeV;
  G4double kinEs[nmax+1];

  G4double cut=10.*keV;
  for (int j=0; j<=nmax; ++j) 
    kinEs[j]=exp(G4double(j)/nmax*log(emax/emin))*emin;
 

  // writing information to root file
  // Z , cut, emin, emax
  G4int currentZ=0;

  TTree tree("info","info");
  tree.Branch("Z",&currentZ,"Z/I");
  tree.Branch("cut",&cut,"cut/D");
  tree.Branch("emin",&emin,"emin/D");
  tree.Branch("emax",&emax,"emax/D");

  // setup test materials
  const G4int nElements=6;
  G4int theZ[]={1,8,29,47,82,92};

  G4Element* els[nElements];
  G4Material* mats[nElements];
  for (G4int i=0; i<nElements; i++) {
    G4Element * el = G4NistManager::Instance()->FindOrBuildElement(theZ[i]);
    std::vector<G4String> names = G4NistManager::Instance()->GetNistMaterialNames(); 
    G4Material* mat =  G4NistManager::Instance()->FindOrBuildMaterial(names[theZ[i]-1]);
    els[i]=el;
    mats[i]=mat;
  }


  // fill plots
  G4double cross[nmax+1], cross1[nmax+1], cross2[nmax+1];
  G4double xDEDX[nmax+1], xDEDX1[nmax+1], xDEDX2[nmax+1];

  for (G4int i=0; i<nElements; i++ ) {
    currentZ=theZ[i];
    tree.Fill();
    G4cout<<" Z="<<theZ[i]<<" ("<<mats[i]->GetName()<<","<<els[i]->GetName()<<")"<<G4endl;
    const G4double* theAtomicNumDensityVector = mats[i]->GetAtomicNumDensityVector();
    G4double dndV = 0.0;
    for (size_t j=0; j<mats[i]->GetNumberOfElements(); j++) 
      dndV += theAtomicNumDensityVector[j]; 

    for (int j=0; j<=nmax; ++j) {

      modelR->SetupForMaterial(0,mats[i],.1*GeV); // WARNING kinetic energy will be overwritten anyway
      modelR->G4VEmModel::SetCurrentElement(els[i]);
      model2->SetupForMaterial(0,mats[i],.1*GeV); // WARNING kinetic energy will be overwritten anyway
      model2->G4VEmModel::SetCurrentElement(els[i]);

      xDEDX[j]= modelR->ComputeDEDXPerVolume(mats[i], G4Electron::Electron(), kinEs[j], cut)/dndV/barn;
      xDEDX1[j] = model1->ComputeDEDXPerVolume(mats[i], G4Electron::Electron(), kinEs[j], cut)/dndV/barn; 
      xDEDX2[j] = model2->ComputeDEDXPerVolume(mats[i], G4Electron::Electron(), kinEs[j], cut)/dndV/barn; 
      cross[j] = modelR->CrossSectionPerVolume(mats[i], G4Electron::Electron(), kinEs[j], cut, kinEs[j])/dndV/barn;
      cross1[j] = model1->CrossSectionPerVolume(mats[i], G4Electron::Electron(), kinEs[j], cut, kinEs[j])/dndV/barn; 
      cross2[j] = model2->CrossSectionPerVolume(mats[i], G4Electron::Electron(), kinEs[j], cut, kinEs[j])/dndV/barn; 
    }
    TGraph gR(nmax,kinEs,cross);
    gR.Draw("Alp");
    gR.SetLineColor(2);
    gR.GetHistogram()->SetXTitle("E_{kin}");
    gR.GetHistogram()->SetYTitle("#sigma(E_{kin})");
    gR.Write("modelR");
    TGraph g1(nmax,kinEs,cross1);
    g1.Draw("lp");
    g1.SetLineColor(4);
    g1.GetHistogram()->SetXTitle("E_{kin}");
    g1.GetHistogram()->SetYTitle("#sigma(E_{kin})");
    g1.Write("model1");
    TGraph g2(nmax,kinEs,cross2);
    g2.Draw("lp");
    g2.SetLineColor(3);
    g2.SetLineStyle(2);
    g2.GetHistogram()->SetXTitle("E_{kin}");
    g2.GetHistogram()->SetYTitle("#sigma(E_{kin})");
    g2.Write("model2");

    TGraph xR(nmax,kinEs,xDEDX);
    xR.Draw("Alp");
    xR.SetLineColor(2);
    xR.GetHistogram()->SetXTitle("E_{kin}");
    xR.GetHistogram()->SetYTitle("dE/dx");
    xR.Write("xDEDXR");
    TGraph x1(nmax,kinEs,xDEDX1);
    x1.Draw("lp");
    x1.SetLineColor(4);
    x1.GetHistogram()->SetXTitle("E_{kin}");
    x1.GetHistogram()->SetYTitle("dE/dx");
    x1.Write("xDEDX1");
    TGraph x2(nmax,kinEs,xDEDX2);
    x2.Draw("lp");
    x2.SetLineColor(3);
    x2.SetLineStyle(2);
    x2.GetHistogram()->SetXTitle("E_{kin}");
    x2.GetHistogram()->SetYTitle("dE/dx");
    x2.Write("xDEDX2");
  }

  // % differential xsection
  G4double ekin = 1.*GeV;
  emin = 10.*MeV;
  G4double egamma = emin;
  G4double dxsec = 0.;

  TTree tree_dxsec("dxsec","dxsec");
  tree_dxsec.Branch("Z",&currentZ,"Z/I");
  tree_dxsec.Branch("ekin",&ekin,"ekin/D");
  tree_dxsec.Branch("egamma",&egamma,"egamma/D");
  tree_dxsec.Branch("dxsec",&dxsec,"dxsec/D");
  std::cout<<" Calc DX"<<std::endl;
  {
    G4Material* mat = mats[2]; // [2] Cu
    currentZ=theZ[2]; // [2]

    const G4int n = 100;

    model2->SetupForMaterial(0,mat,ekin); // WARNING kinetic energy will be overwritten anyway
    model2->G4VEmModel::SetCurrentElement(els[2]); // [2]
    model2->SetCurrentElement(currentZ);

    for (G4int i=0;i<=n;++i) {
      egamma=std::pow(ekin/emin,double(i)/n)*emin;
      dxsec = model2->ComputeDXSectionPerAtom(egamma);
      // need to multiply with prefactor!!
      tree_dxsec.Fill();
    }
  }

  tree_dxsec.Write();
  tree.Write();
}

int main(int argc,char** argv) {
  G4UnitDefinition::BuildUnitsTable();

  TFile tree("eBremRel01.root","RECREATE","eBremRel");

  // make some tests
  CalcCrossSection();

  tree.Write();
  tree.Close();  

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
