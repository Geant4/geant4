// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SetCutTest.cc,v 1.2 1999-12-15 14:51:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "G4Timer.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

//    it tests the G4SetCut method(s) --------------------------------
//    created by L.Urban on 23/05/96 ---------------------------------
//    revised by G.Cosmo on 06/06/96 ---------------------------------
//    added Gamma, MuonMinus & Positron by L.Urban on 12/06/96 -------
//    revised by H.Kurashige on 04/07/96 ---------------------------------

int main()
{
  //-------- set output format -------
    G4cout.setf( G4std::ios::scientific, G4std::ios::floatfield );
    G4cout << G4std::setprecision(4) << G4endl;
  //--- write results to the file setcut.out -----
    G4std::ofstream outFile("setcut.out", G4std::ios::out );
    outFile.setf( G4std::ios:: scientific, G4std::ios::floatfield );
    outFile << G4std::setprecision(4) << G4endl;

  //--------- Material definition ---------

  G4Timer theTimer;
  G4double a, z, ez, density;
  G4String name, symbol;
  G4int nel;

  a=12.01*g/mole;
  density=2.265*g/cm3;
  G4Material* C = new G4Material(name="Carbon   ",z=6.,a,density);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", ez=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", ez=8., a);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air      ", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  a = 28.09*g/mole ;
  density = 2.33*g/cm3 ;
  G4Material* Si = new G4Material(name="Silicon  ", z=14., a, density);

  a = 39.95*g/mole;
  density = 1.782e-3*g/cm3;
  G4Material* Ar = new G4Material(name="Argon gas", z=18., a, density);

  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron     ", z=26., a, density);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead     ", z=82., a, density);

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4Material* apttoMaterial;
  G4String MaterialName;
 
//--------- Particle definition ---------

  G4MuonPlus* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus* theMuonMinus = G4MuonMinus::MuonMinusDefinition();
  G4Electron* theElectron = G4Electron::ElectronDefinition();
  G4Positron* thePositron = G4Positron::PositronDefinition();
  G4Gamma* theGamma = G4Gamma::GammaDefinition();

// ------------------------------------------------------------
    G4double* MuP;
    G4double* MuM;
    G4double* El;
    G4double* Po;
    G4double* Ga;

// tests for ranges (in mm) given in array rgs   ............
  const G4int nrang=28;
  G4double range;
  G4double rgs[nrang] = {
                  1.e-3,2.e-3,5.e-3,1.e-2,2.e-2,5.e-2,1.e-1,2.e-1,5.e-1,
                  1.,2.,5.,1.e1,2.e1,5.e1,1.e2,2.e2,5.e2,1.e3,2.e3,5.e3,
                  1.e4,2.e4,5.e4,1.e5,2.e5,5.e5,1.e6 };

// ***************************************************************
  outFile << "SetCutTest results for mu+,mu-,e+,e-,gamma" << G4endl;
  outFile << "==========================================" << G4endl;
  outFile << G4endl;
  for (G4int irang=0 ; irang<nrang; irang++)
  {
    range=rgs[irang];
    theTimer.Start();
    theMuonPlus->SetCuts(range);
    theMuonMinus->SetCuts(range);
    thePositron->SetCuts(range);
    theElectron->SetCuts(range);
    theGamma->SetCuts(range);
    theTimer.Stop();
    G4cout << "     range/absorption length cut in mm = "
         << range << G4endl;
    outFile << G4endl;
    outFile << G4endl;
    outFile << "          range/absorption length cut in mm = "
         << range << G4endl;
    outFile << G4endl;
    outFile << "          ";
    outFile << "time(SetCut)=" << theTimer.GetUserElapsed() << G4endl;
    outFile << G4endl;
    outFile << "            cuts in kinetic energy (MeV)" << G4endl;

    outFile << G4endl;
    outFile << " material     mu+         mu-         e+          e-        gamma "
            << G4endl;
    outFile << "----------------------------------------------------------------------------"
            << G4endl;
    MuP = (*theMuonPlus).GetCutsInEnergy(); 
    MuM = (*theMuonMinus).GetCutsInEnergy(); 
    Po = (*thePositron).GetCutsInEnergy(); 
    El = (*theElectron).GetCutsInEnergy(); 
    Ga = (*theGamma).GetCutsInEnergy(); 

    for (G4int icut=0; icut<theMaterialTable->length(); icut++)
    {
      apttoMaterial = (*theMaterialTable)[icut];
      MaterialName = apttoMaterial->GetName();
      outFile << MaterialName << "  "
              << MuP[icut] << "  "
              << MuM[icut] << "  "
              << Po[icut] << "  "
              << El[icut] << "  "
              << Ga[icut] << G4endl;
    }
    outFile << "----------------------------------------------------------------------------"
            << G4endl;
  }

  return EXIT_SUCCESS;
}
