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
// $Id: RunAction.cc,v 1.1 2007-06-21 15:18:32 jjacquem Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:detector(det), primary(kin)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //set precision for printing
  G4int prec = G4cout.precision(6);
   
  // get particle 
  G4ParticleDefinition* particle = primary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double charge   = particle->GetPDGCharge();    
  G4double energy   = primary->GetParticleGun()->GetParticleEnergy();
 
  // get material
  G4Material* material = detector->GetMaterial();
  G4String matName     = material->GetName();
  G4double density     = material->GetDensity();
  G4double radl        = material->GetRadlen();  

  G4cout << "\n " << partName << " ("
         << G4BestUnit(energy,"Energy") << ") in " 
	 << material->GetName() << " (density: " 
	 << G4BestUnit(density,"Volumic Mass") << ";   radiation length: "
	 << G4BestUnit(radl,   "Length")       << ")" << G4endl;	 

  // get cuts	 
  GetCuts();
  if (charge != 0.) {
   G4cout << "\n  Range cuts : \t gamma "  
                      << std::setw(8) << G4BestUnit(rangeCut[0],"Length")
          << "\t e- " << std::setw(8) << G4BestUnit(rangeCut[1],"Length");
   G4cout << "\n Energy cuts : \t gamma " 
                      << std::setw(8) << G4BestUnit(energyCut[0],"Energy")
	  << "\t e- " << std::setw(8) << G4BestUnit(energyCut[1],"Energy")
	  << G4endl;
   }
   	  
  // get processList and extract EM processes (but not MultipleScattering)
  G4ProcessVector* plist = particle->GetProcessManager()->GetProcessList();
  G4String procName;
  G4double cut;
  std::vector<G4String> emName;
  std::vector<G4double> enerCut;
  size_t length = plist->size();
  for (size_t j=0; j<length; j++) {
     procName = (*plist)[j]->GetProcessName();
     cut = energyCut[1];
     if ((procName == "eBrem")||(procName == "muBrems")) cut = energyCut[0];
     if (((*plist)[j]->GetProcessType() == fElectromagnetic) &&
         (procName != "msc")) {
       emName.push_back(procName);
       enerCut.push_back(cut);
     }  
  }
  
  // print list of processes
  G4cout << "\n  processes :                ";
  for (size_t j=0; j<emName.size();j++)
    G4cout << "\t" << std::setw(13) << emName[j] << "\t";
  G4cout << "\t" << std::setw(13) <<"total";
      
  //instanciate EmCalculator
  G4EmCalculator emCal;
  //  emCal.SetVerbose(2);
  
  //compute cross section per atom (only for single material)
  if (material->GetNumberOfElements() == 1) {
    G4double Z = material->GetZ();
    G4double A = material->GetA();
     
    std::vector<G4double> sigma0;
    G4double sig, sigtot = 0.;

    for (size_t j=0; j<emName.size();j++) {
      sig = emCal.ComputeCrossSectionPerAtom
                      (energy,particle,emName[j],Z,A,enerCut[j]);
      sigtot += sig; 			     
      sigma0.push_back(sig);	        
    }
    sigma0.push_back(sigtot);

    G4cout << "\n \n  cross section per atom   : ";
    for (size_t j=0; j<sigma0.size();j++) {	     
      G4cout << "\t" << std::setw(13) << G4BestUnit(sigma0[j], "Surface");
    }
    G4cout << G4endl;
  }
    
  //get cross section per volume 
  std::vector<G4double> sigma1;
  std::vector<G4double> sigma2;  
  G4double Sig, Sigtot = 0.;

  for (size_t j=0; j<emName.size();j++) {
    Sig = emCal.GetCrossSectionPerVolume(energy,particle,emName[j],material);
    if (Sig == 0.) Sig = emCal.ComputeCrossSectionPerVolume
		     (energy,particle,emName[j],material,enerCut[j]);
    Sigtot += Sig; 			     
    sigma1.push_back(Sig);
    sigma2.push_back(Sig/density);		        
  }
  sigma1.push_back(Sigtot);
  sigma2.push_back(Sigtot/density);	  
    
  //print cross sections
  G4cout << "\n \n  cross section per volume : ";
  for (size_t j=0; j<sigma1.size();j++) {	     
    G4cout << "\t" << std::setw(13) << sigma1[j]*cm << " cm^-1";
  }
  
  G4cout << "\n  cross section per mass   : ";
  for (size_t j=0; j<sigma2.size();j++) {
    G4cout << "\t" << std::setw(13) << G4BestUnit(sigma2[j], "Surface/Mass");
  }
   
  //print mean free path
  
  G4double lambda;
  
  G4cout << "\n \n  mean free path           : ";
  for (size_t j=0; j<sigma1.size();j++) {
    lambda = DBL_MAX; 
    if (sigma1[j] > 0.) lambda = 1/sigma1[j];
    G4cout << "\t" << std::setw(13) << G4BestUnit( lambda, "Length");
  }
  
  //mean free path (g/cm2)
  G4cout << "\n        (g/cm2)            : ";  
  for (size_t j=0; j<sigma2.size();j++) {
    lambda =  DBL_MAX;
    if (sigma2[j] > 0.) lambda = 1/sigma2[j];   	            
    G4cout << "\t" << std::setw(13) << G4BestUnit( lambda, "Mass/Surface");    
  }
  G4cout << G4endl;
  
  if (charge == 0.) {
    G4cout.precision(prec);
    G4cout << "\n-------------------------------------------------------------\n"
           << G4endl;
    return;
  }
  
  //get stopping power 
  std::vector<G4double> dedx1;
  std::vector<G4double> dedx2;  
  G4double dedx, dedxtot = 0.;

  for (size_t j=0; j<emName.size();j++) {
    dedx = emCal.ComputeDEDX(energy,particle,emName[j],material,enerCut[j]);
    dedx1.push_back(dedx);
    dedx2.push_back(dedx/density);		        
  }
  dedxtot = emCal.GetDEDX(energy,particle,material);
  dedx1.push_back(dedxtot);
  dedx2.push_back(dedxtot/density);	  
    
  //print stopping power
  G4cout << "\n \n  restricted dE/dx         : ";
  for (size_t j=0; j<sigma1.size();j++) {	     
    G4cout << "\t" << std::setw(13) << G4BestUnit(dedx1[j],"Energy/Length");
  }
  
  G4cout << "\n      (MeV/g/cm2)          : ";
  for (size_t j=0; j<sigma2.size();j++) {
  G4cout << "\t" << std::setw(13) << G4BestUnit(dedx2[j],"Energy*Surface/Mass");
  }
  
  //get range from restricted dedx
  G4double range1 = emCal.GetRangeFromRestricteDEDX(energy,particle,material);
  G4double range2 = range1*density;
  
   //get range from full dedx
  G4double Range1 = emCal.GetCSDARange(energy,particle,material);
  G4double Range2 = Range1*density;
  
  //print range
  G4cout << "\n \n  range from restrict dE/dx: " 
         << "\t" << std::setw(8) << G4BestUnit(range1,"Length")
         << " (" << std::setw(8) << G4BestUnit(range2,"Mass/Surface") << ")";
	 
  G4cout << "\n  range from full dE/dx    : " 
         << "\t" << std::setw(8) << G4BestUnit(Range1,"Length")
         << " (" << std::setw(8) << G4BestUnit(Range2,"Mass/Surface") << ")";
	 
  //get transport mean free path (for multiple scattering)
  G4double MSmfp1 = emCal.GetMeanFreePath(energy,particle,"msc",material);
  G4double MSmfp2 = MSmfp1*density;
  
  //print transport mean free path
  G4cout << "\n \n  transport mean free path : " 
         << "\t" << std::setw(8) << G4BestUnit(MSmfp1,"Length")
         << " (" << std::setw(8) << G4BestUnit(MSmfp2,"Mass/Surface") << ")";

  if (particle == G4Electron::Electron()) CriticalEnergy();
  	 
  G4cout << "\n-------------------------------------------------------------\n";
  G4cout << G4endl;
       
 // reset default precision
 G4cout.precision(prec);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ProductionCutsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::GetCuts()
{  
  G4ProductionCutsTable* theCoupleTable =
        G4ProductionCutsTable::GetProductionCutsTable();
	
  size_t numOfCouples = theCoupleTable->GetTableSize();
  const G4MaterialCutsCouple* couple = 0;
  G4int index = 0;
  for (size_t i=0; i<numOfCouples; i++) {
     couple = theCoupleTable->GetMaterialCutsCouple(i);
     if (couple->GetMaterial() == detector->GetMaterial()) {index = i; break;}
  }
  
  rangeCut[0] =
         (*(theCoupleTable->GetRangeCutsVector(idxG4GammaCut)))[index];
  rangeCut[1] =      
         (*(theCoupleTable->GetRangeCutsVector(idxG4ElectronCut)))[index];
  rangeCut[2] =      
         (*(theCoupleTable->GetRangeCutsVector(idxG4PositronCut)))[index]; 

  energyCut[0] =
         (*(theCoupleTable->GetEnergyCutsVector(idxG4GammaCut)))[index];
  energyCut[1] =      
         (*(theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut)))[index];
  energyCut[2] =      
         (*(theCoupleTable->GetEnergyCutsVector(idxG4PositronCut)))[index];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CriticalEnergy()
{
  // compute e- critical energy (Rossi definition) and Moliere radius.
  // Review of Particle Physics - Eur. Phys. J. C3 (1998) page 147
  //
  G4EmCalculator emCal;
    
  const G4Material* material = detector->GetMaterial();
  const G4double radl = material->GetRadlen();
  G4double ekin = 5*MeV;
  G4double deioni;
  G4double err  = 1., errmax = 0.001;
  G4int    iter = 0 , itermax = 10;  
  while (err > errmax && iter < itermax) {
    iter++;          
    deioni  = radl*
              emCal.ComputeDEDX(ekin,G4Electron::Electron(),"eIoni",material);
    err = std::abs(deioni - ekin)/ekin;
    ekin = deioni;
  }
  G4cout << "\n \n  critical energy (Rossi)  : " 
         << "\t" << std::setw(8) << G4BestUnit(ekin,"Energy");
	 
  //Pdg formula (only for single material)
  G4double pdga[2] = { 610*MeV, 710*MeV };
  G4double pdgb[2] = { 1.24, 0.92 };
  G4double EcPdg = 0.;
  
  if (material->GetNumberOfElements() == 1) {
    G4int istat = 0;
    if (material->GetState() == kStateGas) istat = 1;  
    G4double Zeff = material->GetZ() + pdgb[istat];
    EcPdg = pdga[istat]/Zeff;
    G4cout << "\t\t\t (from Pdg formula : " 
           << std::setw(8) << G4BestUnit(EcPdg,"Energy") << ")";    
  }
     
 const G4double Es = 21.2052*MeV;
 G4double rMolier1 = Es/ekin, rMolier2 = rMolier1*radl;
 G4cout << "\n  Moliere radius           : "
        << "\t" << std::setw(8) << rMolier1 << " X0 "   
        << "= " << std::setw(8) << G4BestUnit(rMolier2,"Length");
	
 if (material->GetNumberOfElements() == 1) {
    G4double rMPdg = radl*Es/EcPdg;
    G4cout << "\t (from Pdg formula : " 
           << std::setw(8) << G4BestUnit(rMPdg,"Length") << ")";    
  }	 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
