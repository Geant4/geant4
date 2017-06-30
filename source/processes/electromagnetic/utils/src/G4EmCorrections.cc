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
// $Id: G4EmCorrections.cc 103954 2017-05-04 11:29:22Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class 
//
// File name:     G4EmCorrections
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.01.2005
//
// Modifications:
// 05.05.2005 V.Ivanchenko Fix misprint in Mott term
// 26.11.2005 V.Ivanchenko Fix effective charge for heavy ions using original paper
// 28.04.2006 V.Ivanchenko General cleanup, add finite size corrections
// 13.05.2006 V.Ivanchenko Add corrections for ion stopping
// 08.05.2007 V.Ivanchenko Use G4IonTable for ion mass instead of NistTable to avoid
//                         division by zero
// 29.02.2008 V.Ivanchenko use expantions for log and power function
// 21.04.2008 Updated computations for ions (V.Ivanchenko)
// 20.05.2008 Removed Finite Size correction (V.Ivanchenko)
// 19.04.2012 Fix reproducibility problem (A.Ribon)
//
//
// Class Description:
//
// This class provides calculation of EM corrections to ionisation
//

// -------------------------------------------------------------------
//

#include "G4EmCorrections.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4VEmModel.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4AtomicShells.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"
#include "G4Threading.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4double G4EmCorrections::inveplus = 1.0/CLHEP::eplus;
const G4double G4EmCorrections::ZD[11] = 
    {0., 0., 0., 1.72, 2.09, 2.48, 2.82, 3.16, 3.53, 3.84, 4.15};
const G4double G4EmCorrections::UK[20] = {1.9999, 2.0134, 2.0258, 2.0478, 2.0662,
                           2.0817, 2.0945, 2.0999, 2.1049, 2.1132,
                           2.1197, 2.1246, 2.1280, 2.1292, 2.1301,
                           2.1310, 2.1310, 2.1300, 2.1283, 2.1271};
const G4double G4EmCorrections::VK[20] = {8.3410, 8.3373, 8.3340, 8.3287, 8.3247,
                           8.3219, 8.3201, 8.3194, 8.3191, 8.3188,
                           8.3191, 8.3199, 8.3211, 8.3218, 8.3226,
                           8.3244, 8.3264, 8.3285, 8.3308, 8.3320};
G4double G4EmCorrections::ZK[] = {0.0};
const G4double G4EmCorrections::Eta[29] = {0.005, 0.007, 0.01, 0.015, 0.02,
                            0.03,  0.04,  0.05, 0.06,  0.08,
                            0.1,   0.15,  0.2,  0.3,   0.4,
                            0.5,   0.6,   0.7,  0.8,   1.0,
                            1.2,   1.4,   1.5,  1.7,   2.0, 3.0, 5.0, 7.0, 10.};
G4double G4EmCorrections::CK[20][29];
G4double G4EmCorrections::CL[26][28];
const G4double G4EmCorrections::UL[] = {0.1215, 0.5265, 0.8411, 1.0878, 1.2828,
                           1.4379, 1.5032, 1.5617, 1.6608, 1.7401,
                           1.8036, 1.8543, 1.8756, 1.8945, 1.9262,
                           1.9508, 1.9696, 1.9836, 1.9890, 1.9935,
                           2.0001, 2.0039, 2.0053, 2.0049, 2.0040, 2.0028};  
G4double G4EmCorrections::VL[] = {0.0};

G4LPhysicsFreeVector* G4EmCorrections::BarkasCorr = nullptr;
G4LPhysicsFreeVector* G4EmCorrections::ThetaK = nullptr;
G4LPhysicsFreeVector* G4EmCorrections::ThetaL = nullptr;

G4EmCorrections::G4EmCorrections(G4int verb)
{
  particle   = nullptr;
  curParticle= nullptr;
  material   = nullptr;
  curMaterial= nullptr;
  theElementVector = nullptr;
  atomDensity= nullptr;
  curVector  = nullptr;
  ionLEModel = nullptr;
  ionHEModel = nullptr;

  kinEnergy  = 0.0;
  verbose    = verb;
  massFactor = 1.0;
  eth        = 2.0*CLHEP::MeV;
  nbinCorr   = 20;
  eCorrMin   = 25.*CLHEP::keV;
  eCorrMax   = 250.*CLHEP::MeV;

  ionTable = G4ParticleTable::GetParticleTable()->GetIonTable();
  g4calc = G4Pow::GetInstance();

  nIons = ncouples = numberOfElements = idx = currentZ = 0;
  mass = tau = gamma = bg2 = beta2 = beta = ba2 = tmax = charge = q2 = 0.0;

  // Constants
  alpha2 = CLHEP::fine_structure_const*CLHEP::fine_structure_const;

  // G.S. Khandelwal Nucl. Phys. A116(1968)97 - 111.
  // "Shell corrections for K- and L- electrons

  nK = 20;
  nL = 26;
  nEtaK = 29;
  nEtaL = 28;

  isMaster = false;

  // fill vectors
  if(BarkasCorr == nullptr) { Initialise(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCorrections::~G4EmCorrections()
{
  for(G4int i=0; i<nIons; ++i) {delete stopData[i];}
  if(isMaster) { 
    delete BarkasCorr;
    delete ThetaK;
    delete ThetaL;
    BarkasCorr = ThetaK = ThetaL = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::HighOrderCorrections(const G4ParticleDefinition* p,
                                               const G4Material* mat,
                                               G4double e, G4double)
{
  // . Z^3 Barkas effect in the stopping power of matter for charged particles
  //   J.C Ashley and R.H.Ritchie
  //   Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
  //   and ICRU49 report
  //   valid for kineticEnergy < 0.5 MeV
  //   Other corrections from S.P.Ahlen Rev. Mod. Phys., Vol 52, No1, 1980

  SetupKinematics(p, mat, e);
  if(tau <= 0.0) { return 0.0; }

  G4double Barkas = BarkasCorrection (p, mat, e);
  G4double Bloch  = BlochCorrection (p, mat, e);
  G4double Mott   = MottCorrection (p, mat, e);

  G4double sum = (2.0*(Barkas + Bloch) + Mott);

  if(verbose > 1) {
    G4cout << "EmCorrections: E(MeV)= " << e/MeV << " Barkas= " << Barkas
           << " Bloch= " << Bloch << " Mott= " << Mott 
           << " Sum= " << sum << " q2= " << q2 << G4endl; 
    G4cout << " ShellCorrection: " << ShellCorrection(p, mat, e) 
           << " Kshell= " << KShellCorrection(p, mat, e)
           << " Lshell= " << LShellCorrection(p, mat, e)
           << "   " << mat->GetName() << G4endl;
  }
  sum *= material->GetElectronDensity() * q2 *  twopi_mc2_rcl2 /beta2;
  return sum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::IonBarkasCorrection(const G4ParticleDefinition* p,
                                              const G4Material* mat,
                                              G4double e)
{
  // . Z^3 Barkas effect in the stopping power of matter for charged particles
  //   J.C Ashley and R.H.Ritchie
  //   Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
  //   and ICRU49 report
  //   valid for kineticEnergy < 0.5 MeV

  SetupKinematics(p, mat, e);
  G4double res = 0.0;
  if(tau > 0.0) 
    res = 2.0*BarkasCorrection(p, mat, e)*
      material->GetElectronDensity() * q2 *  twopi_mc2_rcl2 /beta2;
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::ComputeIonCorrections(const G4ParticleDefinition* p,
                                                const G4Material* mat,
                                                G4double e)
{
  // . Z^3 Barkas effect in the stopping power of matter for charged particles
  //   J.C Ashley and R.H.Ritchie
  //   Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
  //   and ICRU49 report
  //   valid for kineticEnergy < 0.5 MeV
  //   Other corrections from S.P.Ahlen Rev. Mod. Phys., Vol 52, No1, 1980
  SetupKinematics(p, mat, e);
  if(tau <= 0.0) { return 0.0; }

  G4double Barkas = BarkasCorrection (p, mat, e);
  G4double Bloch  = BlochCorrection (p, mat, e);
  G4double Mott   = MottCorrection (p, mat, e);

  G4double sum = 2.0*(Barkas*(charge - 1.0)/charge + Bloch) + Mott;

  if(verbose > 1) {
    G4cout << "EmCorrections: E(MeV)= " << e/MeV << " Barkas= " << Barkas
           << " Bloch= " << Bloch << " Mott= " << Mott 
           << " Sum= " << sum << G4endl; 
  }
  sum *= material->GetElectronDensity() * q2 *  twopi_mc2_rcl2 /beta2;

  if(verbose > 1) { G4cout << " Sum= " << sum << G4endl; } 
  return sum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::IonHighOrderCorrections(const G4ParticleDefinition* p,
                                                  const G4MaterialCutsCouple* couple,
                                                  G4double e)
{
  // . Z^3 Barkas effect in the stopping power of matter for charged particles
  //   J.C Ashley and R.H.Ritchie
  //   Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
  //   and ICRU49 report
  //   valid for kineticEnergy < 0.5 MeV
  //   Other corrections from S.P.Ahlen Rev. Mod. Phys., Vol 52, No1, 1980

  G4double sum = 0.0;

  if(ionHEModel) {
    G4int Z = G4lrint(p->GetPDGCharge()*inveplus);
    if(Z >= 100)   Z = 99;
    else if(Z < 1) Z = 1;

    G4double ethscaled = eth*p->GetPDGMass()/proton_mass_c2;
    G4int ionPDG = p->GetPDGEncoding();
    if(thcorr.find(ionPDG)==thcorr.end()) {  // Not found: fill the map
      std::vector<G4double> v;
      for(size_t i=0; i<ncouples; ++i){
        v.push_back(ethscaled*ComputeIonCorrections(p,currmat[i],ethscaled));
      }
      thcorr.insert(std::pair< G4int, std::vector<G4double> >(ionPDG,v)); 
    }

    //G4cout << " map size=" << thcorr.size() << G4endl;
    //for(std::map< G4int, std::vector<G4double> >::iterator 
    //    it = thcorr.begin(); it != thcorr.end(); ++it){
    //  G4cout << "\t map element: first (key)=" << it->first  
    //     << "\t second (vector): vec size=" << (it->second).size() << G4endl;
    //  for(size_t i=0; i<(it->second).size(); ++i){
    // G4cout << "\t \t vec element: [" << i << "]=" << (it->second)[i]
    //<< G4endl; } }

    G4double rest = (thcorr.find(ionPDG)->second)[couple->GetIndex()];

    sum = ComputeIonCorrections(p,couple->GetMaterial(),e) - rest/e;

    if(verbose > 1) { 
      G4cout << " Sum= " << sum << " dSum= " << rest/e << G4endl; 
    } 
  }
  return sum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::Bethe(const G4ParticleDefinition* p,
                                const G4Material* mat, 
                                G4double e)
{
  SetupKinematics(p, mat, e);
  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc2 = eexc*eexc;
  G4double dedx = 0.5*G4Log(2.0*electron_mass_c2*bg2*tmax/eexc2)-beta2;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::SpinCorrection(const G4ParticleDefinition* p,
                                         const G4Material* mat,
                                         G4double e)
{
  SetupKinematics(p, mat, e);
  G4double dedx  = 0.5*tmax/(kinEnergy + mass);
  return 0.5*dedx*dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections:: KShellCorrection(const G4ParticleDefinition* p,
                                            const G4Material* mat, 
                                            G4double e)
{
  SetupKinematics(p, mat, e);
  G4double term = 0.0;
  for (G4int i = 0; i<numberOfElements; ++i) {

    G4double Z = (*theElementVector)[i]->GetZ();
    G4int   iz = G4lrint(Z);
    G4double f = 1.0;
    G4double Z2= (Z-0.3)*(Z-0.3);
    if(1 == iz) {
      f  = 0.5;
      Z2 = 1.0;
    }
    G4double eta = ba2/Z2;
    G4double tet = Z2*(1. + Z2*0.25*alpha2);
    if(11 < iz) { tet = ThetaK->Value(Z); }
    term += f*atomDensity[i]*KShell(tet,eta)/Z;
  }

  term /= material->GetTotNbOfAtomsPerVolume();

  return term;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections:: LShellCorrection(const G4ParticleDefinition* p,
                                            const G4Material* mat, 
                                            G4double e)
{
  SetupKinematics(p, mat, e);
  G4double term = 0.0;
  for (G4int i = 0; i<numberOfElements; ++i) {

    G4double Z = (*theElementVector)[i]->GetZ();
    G4int   iz = G4lrint(Z);
    if(2 < iz) {
      G4double Zeff = Z - ZD[10];
      if(iz < 10) { Zeff = Z - ZD[iz]; }
      G4double Z2= Zeff*Zeff;
      G4double f = 0.125;
      G4double eta = ba2/Z2;
      G4double tet = ThetaL->Value(Z);
      G4int nmax = std::min(4,G4AtomicShells::GetNumberOfShells(iz));
      for(G4int j=1; j<nmax; ++j) {
        G4int ne = G4AtomicShells::GetNumberOfElectrons(iz,j);
        if(15 >= iz) {
          if(3 > j) { tet = 0.25*Z2*(1.0 + 5*Z2*alpha2/16.); }
          else      { tet = 0.25*Z2*(1.0 + Z2*alpha2/16.); }
        }
        //G4cout << " LShell: j= " << j << " ne= " << ne << " e(eV)= " << e/eV
        //       << " ThetaL= " << tet << G4endl;
        term += f*ne*atomDensity[i]*LShell(tet,eta)/Z;
      }
    }
  }

  term /= material->GetTotNbOfAtomsPerVolume();

  return term;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::KShell(G4double tet, G4double eta)
{
  G4double corr = 0.0;

  static const G4double TheK[20] = 
    {0.64, 0.65, 0.66, 0.68, 0.70, 0.72, 0.74, 0.75, 0.76, 0.78,
     0.80, 0.82, 0.84, 0.85, 0.86, 0.88, 0.90, 0.92, 0.94, 0.95};


  G4double x = tet;
  G4int itet = 0;
  G4int ieta = 0;
  if(tet < TheK[0]) { 
    x =  TheK[0]; 
  } else if(tet > TheK[nK-1]) { 
    x =  TheK[nK-1];
    itet = nK-2; 
  } else { 
    itet = Index(x, TheK, nK);
  }
  // assimptotic case
  if(eta >= Eta[nEtaK-1]) {
    corr = 
      (Value(x, TheK[itet], TheK[itet+1], UK[itet], UK[itet+1]) + 
       Value(x, TheK[itet], TheK[itet+1], VK[itet], VK[itet+1])/eta +
       Value(x, TheK[itet], TheK[itet+1], ZK[itet], ZK[itet+1])/(eta*eta))/eta;
  } else {
    G4double y = eta;
    if(eta < Eta[0]) { 
      y =  Eta[0]; 
    } else { 
      ieta = Index(y, Eta, nEtaK);
    }
    corr = Value2(x, y, TheK[itet], TheK[itet+1], Eta[ieta], Eta[ieta+1],
                  CK[itet][ieta], CK[itet+1][ieta], 
                  CK[itet][ieta+1], CK[itet+1][ieta+1]);
    //G4cout << "   x= " <<x<<" y= "<<y<<" tet= " <<TheK[itet]
    //           <<" "<< TheK[itet+1]<<" eta= "<< Eta[ieta]<<" "<< Eta[ieta+1]
    //           <<" CK= " << CK[itet][ieta]<<" "<< CK[itet+1][ieta]
    //           <<" "<< CK[itet][ieta+1]<<" "<< CK[itet+1][ieta+1]<<G4endl;
  }
  //G4cout << "Kshell:  tet= " << tet << " eta= " << eta << "  C= " << corr
  //         << " itet= " << itet << "  ieta= " << ieta <<G4endl;
  return corr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::LShell(G4double tet, G4double eta)
{
  G4double corr = 0.0;

  static const G4double TheL[26] = 
    {0.24, 0.26, 0.28, 0.30, 0.32,   0.34, 0.35, 0.36, 0.38, 0.40,
     0.42, 0.44, 0.45, 0.46, 0.48,   0.50, 0.52, 0.54, 0.55, 0.56,
     0.58, 0.60, 0.62, 0.64, 0.65,   0.66};

  G4double x = tet;
  G4int itet = 0;
  G4int ieta = 0;
  if(tet < TheL[0]) { 
    x =  TheL[0]; 
  } else if(tet > TheL[nL-1]) { 
    x =  TheL[nL-1];
    itet = nL-2; 
  } else { 
    itet = Index(x, TheL, nL);
  }

  // assimptotic case
  if(eta >= Eta[nEtaL-1]) {
    corr = (Value(x, TheL[itet], TheL[itet+1], UL[itet], UL[itet+1])
              + Value(x, TheL[itet], TheL[itet+1], VL[itet], VL[itet+1])/eta
            )/eta;
  } else {
    G4double y = eta;
    if(eta < Eta[0]) { 
      y =  Eta[0]; 
    } else { 
      ieta = Index(y, Eta, nEtaL);
    }
    corr = Value2(x, y, TheL[itet], TheL[itet+1], Eta[ieta], Eta[ieta+1],
                  CL[itet][ieta], CL[itet+1][ieta], 
                  CL[itet][ieta+1], CL[itet+1][ieta+1]);
    //G4cout << "   x= " <<x<<" y= "<<y<<" tet= " <<TheL[itet]
    //           <<" "<< TheL[itet+1]<<" eta= "<< Eta[ieta]<<" "<< Eta[ieta+1]
    //           <<" CL= " << CL[itet][ieta]<<" "<< CL[itet+1][ieta]
    //           <<" "<< CL[itet][ieta+1]<<" "<< CL[itet+1][ieta+1]<<G4endl;
  }
  //G4cout<<"Lshell:  tet= "<<tet<<" eta= "<<eta<<" itet= "<<itet
  //        <<" ieta= "<<ieta<<"  Corr= "<<corr<<G4endl;
  return corr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::ShellCorrectionSTD(const G4ParticleDefinition* p,
                                             const G4Material* mat, 
                                             G4double e)
{
  SetupKinematics(p, mat, e);
  G4double taulim= 8.0*MeV/mass;
  G4double bg2lim= taulim * (taulim+2.0);

  G4double* shellCorrectionVector =
            material->GetIonisation()->GetShellCorrectionVector();
  G4double sh = 0.0;
  G4double x  = 1.0;
  G4double taul  = material->GetIonisation()->GetTaul();

  if ( bg2 >= bg2lim ) {
    for (G4int k=0; k<3; ++k) {
        x *= bg2 ;
        sh += shellCorrectionVector[k]/x;
    }

  } else {
    for (G4int k=0; k<3; ++k) {
        x *= bg2lim ;
        sh += shellCorrectionVector[k]/x;
    }
    sh *= G4Log(tau/taul)/G4Log(taulim/taul);
  }
  sh *= 0.5;
  return sh;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::ShellCorrection(const G4ParticleDefinition* p,
                                          const G4Material* mat,
                                          G4double ekin)
{
  SetupKinematics(p, mat, ekin);

  G4double term = 0.0;
  //G4cout << "### G4EmCorrections::ShellCorrection " << mat->GetName()
  //         << "   " << ekin/MeV << " MeV " << G4endl;
  for (G4int i = 0; i<numberOfElements; ++i) {

    G4double res = 0.0;
    G4double res0 = 0.0;
    G4double Z = (*theElementVector)[i]->GetZ();
    G4int   iz = G4lrint(Z);
    G4double Z2= (Z-0.3)*(Z-0.3);
    G4double f = 1.0;
    if(1 == iz) {
      f  = 0.5;
      Z2 = 1.0;
    }
    G4double eta = ba2/Z2;
    G4double tet = Z2*(1. + Z2*0.25*alpha2);
    if(11 < iz) { tet = ThetaK->Value(Z); }
    res0 = f*KShell(tet,eta);
    res += res0;
    //G4cout << " Z= " << iz << " Shell 0" << " tet= " << tet 
    //       << " eta= " << eta << "  resK= " << res0 << G4endl;
    if(2 < iz) {
      G4double Zeff = Z - ZD[10];
      if(iz < 10) { Zeff = Z - ZD[iz]; }
      Z2= Zeff*Zeff;
      eta = ba2/Z2;
      f = 0.125;
      tet = ThetaL->Value(Z);
      G4int ntot = G4AtomicShells::GetNumberOfShells(iz);
      G4int nmax = std::min(4, ntot);
      G4double norm   = 0.0;
      G4double eshell = 0.0;
      for(G4int j=1; j<nmax; ++j) {
        G4int ne = G4AtomicShells::GetNumberOfElectrons(iz,j);
        if(15 >= iz) {
          if(3 > j) { tet = 0.25*Z2*(1.0 + 5*Z2*alpha2/16.); }
          else      { tet = 0.25*Z2*(1.0 + Z2*alpha2/16.); }
        }
        norm   += ne;
        eshell += tet*ne;
        res0 = f*ne*LShell(tet,eta);
        res += res0;
        //G4cout << " Z= " << iz << " Shell " << j << " Ne= " << ne
        //       << " tet= " << tet << " eta= " << eta 
        //       << "  resL= " << res0 << G4endl;
      }
      if(ntot > nmax) {
        eshell /= norm;

	static const G4double HM[53] = {
	  12.0, 12.0, 12.0, 12.0, 11.9, 11.7, 11.5, 11.2, 10.8, 10.4,
	  10.0,  9.51, 8.97, 8.52, 8.03, 7.46, 6.95, 6.53, 6.18, 5.87, 
	  5.61,  5.39, 5.19, 5.01, 4.86, 4.72, 4.62, 4.53, 4.44, 4.38, 
	  4.32,  4.26, 4.20, 4.15, 4.1,  4.04, 4.00, 3.95, 3.93, 3.91, 
	  3.90,  3.89, 3.89, 3.88, 3.88, 3.88, 3.88, 3.88, 3.89, 3.89, 
	  3.90, 3.92, 3.93 };
	static const G4double HN[31] = {
	  75.5, 61.9, 52.2, 45.1, 39.6, 35.4, 31.9, 29.1, 27.2, 25.8, 
	  24.5, 23.6, 22.7, 22.0, 21.4, 20.9, 20.5, 20.2, 19.9, 19.7, 
	  19.5, 19.3, 19.2, 19.1, 18.4, 18.8, 18.7, 18.6, 18.5, 18.4, 
	  18.2};

        // Add M-shell
        if(28 > iz) {
          res += f*(iz - 10)*LShell(eshell,HM[iz-11]*eta);
        } else if(63 > iz) { 
          res += f*18*LShell(eshell,HM[iz-11]*eta);
        } else {
          res += f*18*LShell(eshell,HM[52]*eta);
        }
        // Add N-shell
        if(32 < iz) {
          if(60 > iz) {
            res += f*(iz - 28)*LShell(eshell,HN[iz-33]*eta);
          } else if(63 > iz) {
            res += 4*LShell(eshell,HN[iz-33]*eta);
          } else {
            res += 4*LShell(eshell,HN[30]*eta);
          }
          // Add O-P-shells
          if(60 < iz) {
            res += f*(iz - 60)*LShell(eshell,150*eta);
          }
        }
      }
    }
    term += res*atomDensity[i]/Z;
  }

  term /= material->GetTotNbOfAtomsPerVolume();
  //G4cout << "#     Shell Correction= " << term << G4endl;
  return term;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::DensityCorrection(const G4ParticleDefinition* p,
                                            const G4Material* mat,
                                            G4double e)
{
  SetupKinematics(p, mat, e);

  G4double cden  = material->GetIonisation()->GetCdensity();
  G4double mden  = material->GetIonisation()->GetMdensity();
  G4double aden  = material->GetIonisation()->GetAdensity();
  G4double x0den = material->GetIonisation()->GetX0density();
  G4double x1den = material->GetIonisation()->GetX1density();

  G4double dedx = 0.0;

  // density correction
  static const G4double twoln10 = 2.0*G4Log(10.0);
  G4double x = G4Log(bg2)/twoln10;
  if ( x >= x0den ) {
    dedx = twoln10*x - cden ;
    if ( x < x1den ) dedx += aden*G4Exp(G4Log(x1den-x)*mden) ;
  }

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::BarkasCorrection(const G4ParticleDefinition* p,
                                           const G4Material* mat, 
                                           G4double e)
{
  // . Z^3 Barkas effect in the stopping power of matter for charged particles
  //   J.C Ashley and R.H.Ritchie
  //   Physical review B Vol.5 No.7 1 April 1972 pp. 2393-2397
  //   valid for kineticEnergy > 0.5 MeV

  SetupKinematics(p, mat, e);
  G4double BarkasTerm = 0.0;

  for (G4int i = 0; i<numberOfElements; ++i) {

    G4double Z = (*theElementVector)[i]->GetZ();
    G4int iz = G4lrint(Z);
    if(iz == 47) {
      BarkasTerm += atomDensity[i]*0.006812*G4Exp(-G4Log(beta)*0.9);
    } else if(iz >= 64) {
      BarkasTerm += atomDensity[i]*0.002833*G4Exp(-G4Log(beta)*1.2);
    } else {    

      G4double X = ba2 / Z;
      G4double b = 1.3;
      if(1 == iz) {
        if(material->GetName() == "G4_lH2") { b = 0.6; }
        else                                { b = 1.8; }
      }
      else if(2 == iz)  { b = 0.6; }
      else if(10 >= iz) { b = 1.8; }
      else if(17 >= iz) { b = 1.4; }
      else if(18 == iz) { b = 1.8; }
      else if(25 >= iz) { b = 1.4; }
      else if(50 >= iz) { b = 1.35;}

      G4double W = b/std::sqrt(X);

      G4double val = BarkasCorr->Value(W);
      if(W > BarkasCorr->Energy(46)) { 
        val *= BarkasCorr->Energy(46)/W; 
      } 
      //    G4cout << "i= " << i << " b= " << b << " W= " << W 
      // << " Z= " << Z << " X= " << X << " val= " << val<< G4endl;
      BarkasTerm += val*atomDensity[i] / (std::sqrt(Z*X)*X);
    }
  }

  BarkasTerm *= 1.29*charge/material->GetTotNbOfAtomsPerVolume();

  return BarkasTerm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::BlochCorrection(const G4ParticleDefinition* p,
                                          const G4Material* mat,
                                          G4double e)
{
  SetupKinematics(p, mat, e);

  G4double y2 = q2/ba2;

  G4double term = 1.0/(1.0 + y2);
  G4double del;
  G4double j = 1.0;
  do {
    j += 1.0;
    del = 1.0/(j* (j*j + y2));
    term += del;
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (del > 0.01*term);

  G4double res = -y2*term;
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCorrections::MottCorrection(const G4ParticleDefinition* p,
                                         const G4Material* mat, 
                                         G4double e)
{
  SetupKinematics(p, mat, e);
  G4double mterm = CLHEP::pi*fine_structure_const*beta*charge;
  return mterm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EmCorrections::EffectiveChargeCorrection(const G4ParticleDefinition* p,
                                           const G4Material* mat,
                                           G4double ekin)
{
  G4double factor = 1.0;
  if(p->GetPDGCharge() <= 2.5*CLHEP::eplus || nIons <= 0) { return factor; }
  
  if(verbose > 1) {
    G4cout << "EffectiveChargeCorrection: " << p->GetParticleName() 
           << " in " << mat->GetName()
           << " ekin(MeV)= " << ekin/MeV << G4endl;
  }

  if(p != curParticle || mat != curMaterial) {
    curParticle = p;
    curMaterial = mat;
    curVector   = nullptr;
    currentZ = p->GetAtomicNumber();
    if(verbose > 1) {
      G4cout << "G4EmCorrections::EffectiveChargeCorrection: Zion= " 
             << currentZ << " Aion= " << p->GetPDGMass()/amu_c2 << G4endl;
    }
    massFactor = proton_mass_c2/p->GetPDGMass();
    idx = -1;

    for(G4int i=0; i<nIons; ++i) {
      if(materialList[i] == mat && currentZ == Zion[i]) {
        idx = i;
        break;
      }
    }
    //G4cout << " idx= " << idx << " dz= " << G4endl;
    if(idx >= 0) {
      if(!ionList[idx]) { BuildCorrectionVector(); } 
      if(ionList[idx])  { curVector = stopData[idx]; }
    } else { return factor; }
  }
  if(curVector) {
    factor = curVector->Value(ekin*massFactor);
    if(verbose > 1) {
      G4cout << "E= " << ekin << " factor= " << factor << " massfactor= " 
             << massFactor << G4endl;
    }
  }
  return factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCorrections::AddStoppingData(G4int Z, G4int A,
                                      const G4String& mname,
                                      G4PhysicsVector* dVector)
{
  G4int i = 0;
  for(; i<nIons; ++i) {
    if(Z == Zion[i] && A == Aion[i] && mname == materialName[i]) break;
  }
  if(i == nIons) {
    Zion.push_back(Z);
    Aion.push_back(A);
    materialName.push_back(mname);
    materialList.push_back(nullptr);
    ionList.push_back(nullptr);
    stopData.push_back(dVector);
    nIons++;
    if(verbose > 1) {
      G4cout << "AddStoppingData Z= " << Z << " A= " << A << " " << mname
             << "  idx= " << i << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCorrections::BuildCorrectionVector()
{
  if(!ionLEModel || !ionHEModel) {
    return;
  }

  const G4ParticleDefinition* ion = curParticle;
  G4int Z = Zion[idx];
  if(currentZ != Z) {
    ion = ionTable->GetIon(Z, Aion[idx], 0);
  }
  //G4cout << "BuildCorrectionVector: idx= " << idx << " Z= " << Z 
  //       << " curZ= " << currentZ << G4endl;

  G4double A = G4double(ion->GetBaryonNumber());
  G4PhysicsVector* v = stopData[idx];
    
  const G4ParticleDefinition* p = G4GenericIon::GenericIon();
  G4double massRatio = proton_mass_c2/ion->GetPDGMass();

  if(verbose > 1) {
    G4cout << "BuildCorrectionVector: Stopping for "
           << curParticle->GetParticleName() << " in " 
           << materialName[idx] << " Ion Z= " << Z << " A= " << A
           << " massRatio= " << massRatio << G4endl;
  }

  G4PhysicsLogVector* vv = 
    new G4PhysicsLogVector(eCorrMin,eCorrMax,nbinCorr);
  vv->SetSpline(true);
  G4double e, eion, dedx, dedx1;
  G4double eth0 = v->Energy(0);
  G4double escal = eth/massRatio;
  G4double qe = 
    effCharge.EffectiveChargeSquareRatio(ion, curMaterial, escal); 
  G4double dedxt = 
    ionLEModel->ComputeDEDXPerVolume(curMaterial, p, eth, eth)*qe;
  G4double dedx1t = 
    ionHEModel->ComputeDEDXPerVolume(curMaterial, p, eth, eth)*qe 
    + ComputeIonCorrections(curParticle, curMaterial, escal);
  G4double rest = escal*(dedxt - dedx1t);
  //G4cout << "Escal(MeV)= "<<escal<<" dedxt0= " <<dedxt 
  //       << " dedxt1= " << dedx1t << G4endl;   

  for(G4int i=0; i<=nbinCorr; ++i) {
    e = vv->Energy(i);
    escal = e/massRatio;
    eion  = escal/A;
    if(eion <= eth0) {
      dedx = v->Value(eth0)*std::sqrt(eion/eth0);
    } else {
      dedx = v->Value(eion);
    }
    qe = effCharge.EffectiveChargeSquareRatio(curParticle,curMaterial,escal); 
    if(e <= eth) {
      dedx1 = ionLEModel->ComputeDEDXPerVolume(curMaterial, p, e, e)*qe;
    } else {
      dedx1 = ionHEModel->ComputeDEDXPerVolume(curMaterial, p, e, e)*qe +
        ComputeIonCorrections(curParticle, curMaterial, escal) + rest/escal;
    }
    vv->PutValue(i, dedx/dedx1);
    if(verbose > 1) {
      G4cout << "  E(meV)= " << e/MeV << "   Correction= " << dedx/dedx1
             << "   "  << dedx << " " << dedx1 
             << "  massF= " << massFactor << G4endl;
    }
  }
  delete v;
  ionList[idx]  = ion;
  stopData[idx] = vv;
  if(verbose > 1) { G4cout << "End data set " << G4endl; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCorrections::InitialiseForNewRun()
{
  G4ProductionCutsTable* tb = G4ProductionCutsTable::GetProductionCutsTable();
  ncouples = tb->GetTableSize();
  if(currmat.size() != ncouples) {
    currmat.resize(ncouples);
    for(std::map< G4int, std::vector<G4double> >::iterator it = 
        thcorr.begin(); it != thcorr.end(); ++it){
      (it->second).clear();
    }
    thcorr.clear();
    for(size_t i=0; i<ncouples; ++i) {
      currmat[i] = tb->GetMaterialCutsCouple(i)->GetMaterial();
      G4String nam = currmat[i]->GetName();
      for(G4int j=0; j<nIons; ++j) {
        if(nam == materialName[j]) { materialList[j] = currmat[i]; }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCorrections::Initialise()
{
  if(G4Threading::IsMasterThread()) { isMaster = true; }

  // Z^3 Barkas effect in the stopping power of matter for charged particles
  // J.C Ashley and R.H.Ritchie
  // Physical review B Vol.5 No.7 1 April 1972 pagg. 2393-2397
  G4int i, j;
  static const G4double fTable[47][2] = {
   { 0.02, 21.5},
   { 0.03, 20.0},
   { 0.04, 18.0},
   { 0.05, 15.6},
   { 0.06, 15.0},
   { 0.07, 14.0},
   { 0.08, 13.5},
   { 0.09, 13.},
   { 0.1,  12.2},
   { 0.2,  9.25},
   { 0.3,  7.0},
   { 0.4,  6.0},
   { 0.5,  4.5},
   { 0.6,  3.5},
   { 0.7,  3.0},
   { 0.8,  2.5},
   { 0.9,  2.0},
   { 1.0,  1.7},
   { 1.2,  1.2},
   { 1.3,  1.0},
   { 1.4,  0.86},
   { 1.5,  0.7},
   { 1.6,  0.61},
   { 1.7,  0.52},
   { 1.8,  0.5},
   { 1.9,  0.43},
   { 2.0,  0.42},
   { 2.1,  0.3},
   { 2.4,  0.2},
   { 3.0,  0.13},
   { 3.08, 0.1},
   { 3.1,  0.09},
   { 3.3,  0.08},
   { 3.5,  0.07},
   { 3.8,  0.06},
   { 4.0,  0.051},
   { 4.1,  0.04},
   { 4.8,  0.03},
   { 5.0,  0.024},
   { 5.1,  0.02},
   { 6.0,  0.013},
   { 6.5,  0.01},
   { 7.0,  0.009},
   { 7.1,  0.008},
   { 8.0,  0.006},
   { 9.0,  0.0032},
   { 10.0, 0.0025} };

  BarkasCorr = new G4LPhysicsFreeVector(47, 0.02, 10.);
  for(i=0; i<47; ++i) { BarkasCorr->PutValues(i, fTable[i][0], fTable[i][1]); }
  BarkasCorr->SetSpline(true);

  static const G4double SK[20] = {1.9477, 1.9232, 1.8996, 1.8550, 1.8137,
                           1.7754, 1.7396, 1.7223, 1.7063, 1.6752,
                           1.6461, 1.6189, 1.5933, 1.5811, 1.5693,
                           1.5467, 1.5254, 1.5053, 1.4863, 1.4772};
  static const G4double TK[20] = {2.5222, 2.5125, 2.5026, 2.4821, 2.4608,
                           2.4388, 2.4163, 2.4044, 2.3933, 2.3701,
                           2.3466, 2.3229, 2.2992, 2.2872, 2.2753,
                           2.2515, 2.2277, 2.2040, 2.1804, 2.1686};

  static const G4double SL[26] = {15.3343, 13.9389, 12.7909, 11.8343, 11.0283,
                           10.3424, 10.0371,  9.7537,  9.2443,  8.8005,
                           8.4114,  8.0683,   7.9117, 7.7641, 7.4931,
                           7.2506,  7.0327,   6.8362, 6.7452, 6.6584,
                           6.4969,  6.3498,   6.2154, 6.0923, 6.0345, 5.9792};
  static const G4double TL[26] = {35.0669, 33.4344, 32.0073, 30.7466, 29.6226,
                           28.6128, 28.1449, 27.6991, 26.8674, 26.1061,
                           25.4058, 24.7587, 24.4531, 24.1583, 23.5992,
                           23.0771, 22.5880, 22.1285, 21.9090, 21.6958,
                           21.2872, 20.9006, 20.5341, 20.1859, 20.0183, 19.8546};

  static const G4double bk1[29][11] = { 
  {0.005, 1.34782E-8, 1.46132E-8, 1.72179E-8, 2.03521E-8, 2.41370E-8, 2.87247E-8, 3.13778E-8, 3.43072E-8, 4.11274E-8, 4.94946E-8}, 
  {0.007, 6.87555E-8, 7.44373E-8, 8.74397E-8, 1.03022E-7, 1.21760E-7, 1.44370E-7, 1.57398E-7, 1.71747E-7, 2.05023E-7, 2.45620E-7}, 
  {0.01, 3.78413E-7, 4.08831E-7, 4.78154E-7, 5.60760E-7, 6.59478E-7, 7.77847E-7, 8.45709E-7, 9.20187E-7, 1.09192E-6, 1.29981E-6}, 
  {0.015, 2.53200E-6, 2.72664E-6, 3.16738E-6, 3.68791E-6, 4.30423E-6, 5.03578E-6, 5.45200E-6, 5.90633E-6, 6.94501E-6, 8.18757E-6}, 
  {0.02, 9.43891E-6, 1.01339E-5, 1.16984E-5, 1.35313E-5, 1.56829E-5, 1.82138E-5, 1.96439E-5, 2.11973E-5, 2.47216E-5, 2.88935E-5}, 
  {0.03, 5.67088E-5, 6.05576E-5, 6.91311E-5, 7.90324E-5, 9.04832E-5, 1.03744E-4, 1.11147E-4, 1.19122E-4, 1.36980E-4, 1.57744E-4}, 
  {0.04, 1.91576E-4, 2.03626E-4, 2.30230E-4, 2.60584E-4, 2.95248E-4, 3.34870E-4, 3.56771E-4, 3.80200E-4, 4.32104E-4, 4.91584E-4}, 
  {0.05, 4.74226E-4, 5.02006E-4, 5.62872E-4, 6.31597E-4, 7.09244E-4, 7.97023E-4, 8.45134E-4, 8.96410E-4, 1.00867E-3, 1.13590E-3}, 
  {0.06, 9.67285E-4, 1.02029E-3, 1.13566E-3, 1.26476E-3, 1.46928E-3, 1.57113E-3, 1.65921E-3, 1.75244E-3, 1.95562E-3, 2.18336E-3}, 
  {0.08, 2.81537E-3, 2.95200E-3, 3.24599E-3, 3.57027E-3, 3.92793E-3, 4.32246E-3, 4.53473E-3, 4.75768E-3, 5.23785E-3, 5.76765E-3}, 
  {0.1, 6.14216E-3, 6.40864E-3, 6.97750E-3, 7.59781E-3, 8.27424E-3, 9.01187E-3, 9.40534E-3, 9.81623E-3, 1.06934E-2, 1.16498E-2}, 
  {0.15, 2.23096E-2, 2.30790E-2, 2.46978E-2, 2.64300E-2, 2.82832E-2, 3.02661E-2, 3.13090E-2, 3.23878E-2, 3.46580E-2, 3.70873E-2}, 
  {0.2, 5.04022E-2, 5.18492E-2, 5.48682E-2, 5.80617E-2, 6.14403E-2, 6.50152E-2, 6.68798E-2, 6.87981E-2, 7.28020E-2, 7.70407E-2}, 
  {0.3, 1.38018E-1, 1.41026E-1, 1.47244E-1, 1.53743E-1, 1.60536E-1, 1.67641E-1, 1.71315E-1, 1.75074E-1, 1.82852E-1, 1.90997E-1}, 
  {0.4, 2.56001E-1, 2.60576E-1, 2.69992E-1, 2.79776E-1, 2.89946E-1, 3.00525E-1, 3.05974E-1, 3.11533E-1, 3.22994E-1, 3.34935E-1}, 
  {0.5, 3.92181E-1, 3.98213E-1, 4.10597E-1, 4.23427E-1, 4.36726E-1, 4.50519E-1, 4.57610E-1, 4.64834E-1, 4.79700E-1, 4.95148E-1}, 
  {0.6, 5.37913E-1, 5.45268E-1, 5.60350E-1, 5.75948E-1, 5.92092E-1, 6.08811E-1, 6.17396E-1, 6.26136E-1, 6.44104E-1, 6.62752E-1}, 
  {0.7, 6.87470E-1, 6.96021E-1, 7.13543E-1, 7.31650E-1, 7.50373E-1, 7.69748E-1, 7.79591E-1, 7.89811E-1, 8.10602E-1, 8.32167E-1}, 
  {0.8, 8.37159E-1, 8.46790E-1, 8.66525E-1, 8.86910E-1, 9.07979E-1, 9.29772E-1, 9.40953E-1, 9.52331E-1, 9.75701E-1, 9.99934E-1}, 
  {1.0, 1.12850, 1.14002, 1.16362, 1.18799, 1.21317, 1.23921, 1.25257, 1.26616, 1.29408, 1.32303}, 
  {1.2, 1.40232, 1.41545, 1.44232, 1.47007, 1.49875, 1.52842, 1.54364, 1.55913, 1.59095, 1.62396}, 
  {1.4, 1.65584, 1.67034, 1.70004, 1.73072, 1.76244, 1.79526, 1.81210, 1.82925, 1.86448, 1.90104}, 
  {1.5, 1.77496, 1.79009, 1.82107, 1.85307, 1.88617, 1.92042, 1.93800, 1.95590, 1.99269, 2.03087}, 
  {1.7, 1.99863, 2.01490, 2.04822, 2.08265, 2.11827, 2.15555, 2.17409, 2.19337, 2.23302, 2.27419}, 
  {2.0, 2.29325, 2.31100, 2.34738, 2.38501, 2.42395, 2.46429, 2.48401, 2.50612, 2.54955, 2.59468}, 
  {3.0, 3.08887, 3.11036, 3.15445, 3.20013, 3.24748, 3.29664, 3.32192, 3.34770, 3.40081, 3.45611}, 
  {5.0, 4.07599, 4.10219, 4.15606, 4.21199, 4.27010, 4.33056, 4.36172, 4.39353, 4.45918, 4.52772}, 
  {7.0, 4.69647, 4.72577, 4.78608, 4.84877, 4.91402, 4.98200, 5.01707, 5.05290, 5.12695, 5.20436}, 
  {10.0, 5.32590, 5.35848, 5.42560, 5.49547, 5.56830, 5.64429, 5.68353, 5.72366, 5.80666, 5.89359}
  }; 

  static const G4double bk2[29][11] = { 
  {0.005, 5.98040E-8, 7.25636E-8, 8.00602E-8, 8.84294E-8, 1.08253E-7, 1.33148E-7, 1.64573E-7, 2.04459E-7, 2.28346E-7, 2.55370E-7}, 
  {0.007, 2.95345E-7, 3.56497E-7, 3.92247E-7, 4.32017E-7, 5.25688E-7, 6.42391E-7, 7.88464E-7, 9.72171E-7, 1.08140E-6, 1.20435E-6}, 
  {0.01, 1.55232E-6, 1.86011E-6, 2.03881E-6, 2.23662E-6, 2.69889E-6, 3.26860E-6, 3.26860E-6, 4.84882E-6, 5.36428E-6, 5.94048E-6}, 
  {0.015, 9.67802E-6, 1.14707E-5, 1.25008E-5, 1.36329E-5, 1.62480E-5, 1.94200E-5, 2.32783E-5, 2.79850E-5, 3.07181E-5, 3.37432E-5}, 
  {0.02, 3.38425E-5, 3.97259E-5, 4.30763E-5, 4.67351E-5, 5.51033E-5, 6.51154E-5, 7.71154E-5, 9.15431E-5, 9.98212E-5, 1.08909E-4}, 
  {0.03, 1.81920E-4, 2.10106E-4, 2.25918E-4, 2.43007E-4, 2.81460E-4, 3.26458E-4, 3.79175E-4, 4.41006E-4, 4.75845E-4, 5.13606E-4}, 
  {0.04, 5.59802E-4, 6.38103E-4, 6.81511E-4, 7.28042E-4, 8.31425E-4, 9.50341E-4, 1.08721E-3, 1.24485E-3, 1.33245E-3, 1.42650E-3}, 
  {0.05, 1.28002E-3, 1.44336E-3, 1.53305E-3, 1.62855E-3, 1.83861E-3, 2.07396E-3, 2.34750E-3, 2.65469E-3, 2.82358E-3, 3.00358E-3}, 
  {0.06, 2.42872E-3, 2.72510E-3, 2.88111E-3, 3.04636E-3, 3.40681E-3, 3.81132E-3, 4.26536E-3, 4.77507E-3, 5.05294E-3, 5.34739E-3}, 
  {0.08, 6.35222E-3, 6.99730E-3, 7.34446E-3, 7.70916E-3, 8.49478E-3, 9.36187E-3, 1.03189E-2, 1.13754E-2, 1.19441E-2, 1.25417E-2}, 
  {0.1, 1.26929E-2, 1.38803E-2, 1.44371E-2, 1.50707E-2, 1.64235E-2, 1.78989E-2, 1.95083E-2, 2.12640E-2, 2.22009E-2, 2.31795E-2}, 
  {0.15, 3.96872E-2, 4.24699E-2, 4.39340E-2, 4.54488E-2, 4.86383E-2, 5.20542E-2, 5.57135E-2, 5.96350E-2, 6.17003E-2, 6.38389E-2}, 
  {0.2, 8.15290E-2, 8.62830E-2, 8.87650E-2, 9.13200E-2, 9.66589E-2, 1.02320E-1, 1.08326E-1, 1.14701E-1, 1.18035E-1, 1.21472E-1}, 
  {0.3, 1.99528E-1, 2.08471E-1, 2.13103E-1, 2.17848E-1, 2.27689E-1, 2.38022E-1, 2.48882E-1, 2.60304E-1, 2.66239E-1, 2.72329E-1}, 
  {0.4, 3.47383E-1, 3.60369E-1, 3.67073E-1, 3.73925E-1, 3.88089E-1, 4.02900E-1, 4.18404E-1, 4.34647E-1, 4.43063E-1, 4.51685E-1}, 
  {0.5, 5.11214E-1, 5.27935E-1, 5.36554E-1, 5.45354E-1, 5.63515E-1, 5.82470E-1, 6.02273E-1, 6.22986E-1, 6.33705E-1, 6.44677E-1}, 
  {0.6, 6.82122E-1, 7.02260E-1, 7.12631E-1, 7.23214E-1, 7.45041E-1, 7.67800E-1, 7.91559E-1, 8.16391E-1, 8.29235E-1, 8.42380E-1}, 
  {0.7, 8.54544E-1, 8.77814E-1, 8.89791E-1, 9.02008E-1, 9.27198E-1, 9.53454E-1, 9.80856E-1, 1.00949, 1.02430, 1.03945}, 
  {0.8, 1.02508, 1.05121, 1.06466, 1.07838, 1.10667, 1.13615, 1.16692, 1.19907, 1.21570, 1.23272}, 
  {1.0, 1.35307, 1.38429, 1.40036, 1.41676, 1.45057, 1.48582, 1.52263, 1.56111, 1.58102, 1.60140}, 
  {1.2, 1.65823, 1.69385, 1.71220, 1.73092, 1.76954, 1.80983, 1.85192, 1.89596, 1.91876, 1.94211}, 
  {1.4, 1.93902, 1.97852, 1.99887, 2.01964, 2.06251, 2.10727, 2.15406, 2.20304, 2.22842, 2.25442}, 
  {1.5, 2.07055, 2.11182, 2.13309, 2.15480, 2.19963, 2.24644, 2.29539, 2.34666, 2.37323, 2.40045}, 
  {1.7, 2.31700, 2.36154, 2.38451, 2.40798, 2.45641, 2.50703, 2.56000, 2.61552, 2.64430, 2.67381}, 
  {2.0, 2.64162, 2.69053, 2.71576, 2.74154, 2.79481, 2.85052, 2.90887, 2.97009, 3.00185, 3.03442}, 
  {3.0, 3.51376, 3.57394, 3.60503, 3.63684, 3.70268, 3.77170, 3.84418, 3.92040, 3.96003, 4.00073}, 
  {5.0, 4.59935, 4.67433, 4.71316, 4.75293, 4.83543, 4.92217, 5.01353, 5.10992, 5.16014, 5.21181}, 
  {7.0, 5.28542, 5.37040, 5.41445, 5.45962, 5.55344, 5.65226, 5.79496, 5.90517, 5.96269, 6.02191}, 
  {10.0, 5.98474, 6.08046, 6.13015, 6.18112, 6.28715, 6.39903, 6.51728, 6.64249, 6.70792, 6.77535}
  }; 

  static const G4double bls1[28][10] = { 
  {0.005, 2.4111E-4, 2.5612E-4, 2.7202E-4, 3.0658E-4, 3.4511E-4, 3.8795E-4, 4.3542E-4, 4.6100E-4, 4.8786E-4}, 
  {0.007, 6.3947E-4, 6.7058E-4, 7.0295E-4, 7.7167E-4, 8.4592E-4, 9.2605E-4, 1.0125E-3, 1.0583E-3, 1.1058E-3}, 
  {0.01, 1.5469E-3, 1.6036E-3, 1.6622E-3, 1.7856E-3, 1.9181E-3, 2.1615E-3, 2.3178E-3, 2.4019E-3, 2.4904E-3}, 
  {0.015, 3.7221E-3, 3.8404E-3, 3.9650E-3, 4.2354E-3, 4.5396E-3, 4.8853E-3, 5.2820E-3, 5.5031E-3, 5.7414E-3}, 
  {0.02, 6.9449E-3, 7.1910E-3, 7.4535E-3, 8.0336E-3, 8.6984E-3, 9.4638E-3, 1.0348E-2, 1.0841E-2, 1.1372E-2}, 
  {0.03, 1.7411E-2, 1.8180E-2, 1.8997E-2, 2.0784E-2, 2.2797E-2, 2.5066E-2, 2.7622E-2, 2.9020E-2, 3.0503E-2}, 
  {0.04, 3.8474E-2, 4.0056E-2, 4.1718E-2, 4.5300E-2, 4.9254E-2, 5.3619E-2, 5.8436E-2, 6.1028E-2, 6.3752E-2}, 
  {0.05, 6.7371E-2, 6.9928E-2, 7.2596E-2, 7.8282E-2, 8.4470E-2, 9.1206E-2, 9.8538E-2, 1.0244E-1, 1.0652E-1}, 
  {0.06, 1.0418E-1, 1.0778E-1, 1.1152E-1, 1.1943E-1, 1.2796E-1, 1.3715E-1, 1.4706E-1, 1.5231E-1, 1.5776E-1}, 
  {0.08, 1.9647E-1, 2.0217E-1, 2.0805E-1, 2.2038E-1, 2.3351E-1, 2.4751E-1, 2.6244E-1, 2.7027E-1, 2.7837E-1}, 
  {0.1, 3.0594E-1, 3.1361E-1, 3.2150E-1, 3.3796E-1, 3.5537E-1, 3.7381E-1, 3.9336E-1, 4.0357E-1, 4.1410E-1}, 
  {0.15, 6.1411E-1, 6.2597E-1, 6.3811E-1, 6.6330E-1, 6.8974E-1, 7.1753E-1, 7.4678E-1, 7.6199E-1, 7.7761E-1}, 
  {0.2, 9.3100E-1, 9.5614E-1, 9.7162E-1, 1.0037, 1.0372, 1.0723, 1.1092, 1.1284, 1.1480}, 
  {0.3, 1.5172, 1.5372, 1.5576, 1.5998, 1.6438, 1.6899, 1.7382, 1.7632, 1.7889}, 
  {0.4, 2.0173, 2.0408, 2.0647, 2.1142, 2.1659, 2.2199, 2.2765, 2.3059, 2.3360}, 
  {0.5, 2.3932, 2.4193, 2.4460, 2.5011, 2.5587, 2.6190, 2.6821, 2.7148, 2.7484}, 
  {0.6, 2.7091, 2.7374, 2.7663, 2.8260, 2.8884, 2.9538, 3.0222, 3.0577, 3.0941}, 
  {0.7, 2.9742, 3.0044, 3.0352, 3.0988, 3.1652, 3.2349, 3.3079, 3.3457, 3.3845}, 
  {0.8, 3.2222, 3.2539, 3.2863, 3.3532, 3.4232, 3.4965, 3.5734, 3.6133, 3.6542}, 
  {1.0, 3.6690, 3.7033, 3.7384, 3.8108, 3.8866, 3.9661, 4.0495, 4.0928, 4.1371}, 
  {1.2, 3.9819, 4.0183, 4.0555, 4.1324, 4.2130, 4.2974, 4.3861, 4.4321, 4.4794}, 
  {1.4, 4.2745, 4.3127, 4.3517, 4.4324, 4.5170, 4.6056, 4.6988, 4.7471, 4.7968}, 
  {1.5, 4.4047, 4.4436, 4.4834, 4.5658, 4.6521, 4.7426, 4.8378, 4.8872, 4.9379}, 
  {1.7, 4.6383, 4.6787, 4.7200, 4.8054, 4.8949, 4.9888, 5.0876, 5.1388, 5.1915}, 
  {2.0, 4.9369, 4.9791, 5.0223, 5.1116, 5.2053, 5.3036, 5.4070, 5.4607, 5.5159}, 
  {3.0, 5.6514, 5.6981, 5.7459, 5.8450, 5.9489, 6.0581, 6.1730, 6.2328, 6.2943}, 
  {5.0, 6.4665, 6.5189, 6.5724, 6.6835, 6.8003, 6.9231, 7.0525, 7.1199, 7.1892}, 
  {7.0, 6.8634, 6.9194, 6.9767, 7.0957, 7.2208, 7.3526, 7.4915, 7.5639, 7.6384}
  };
 
  static const G4double bls2[28][10] = { 
  {0.005, 5.4561E-4, 6.0905E-4, 6.7863E-4, 7.5494E-4, 7.9587E-4, 8.3883E-4, 9.3160E-4, 1.0352E-3, 1.1529E-3}, 
  {0.007, 1.2068E-3, 1.3170E-3, 1.4377E-3, 1.5719E-3, 1.6451E-3, 1.7231E-3, 1.8969E-3, 2.1009E-3, 2.3459E-3}, 
  {0.01, 2.6832E-3, 2.9017E-3, 3.1534E-3, 3.4479E-3, 3.6149E-3, 3.7976E-3, 4.2187E-3, 4.7320E-3, 5.3636E-3}, 
  {0.015, 6.2775E-3, 6.9077E-3, 7.6525E-3, 8.5370E-2, 9.0407E-3, 9.5910E-3, 1.0850E-2, 1.2358E-2, 1.4165E-2}, 
  {0.02, 1.2561E-2, 1.3943E-2, 1.5553E-2, 1.7327E-2, 1.8478E-2, 1.9612E-2, 2.2160E-2, 2.5130E-2, 2.8594E-2}, 
  {0.03, 3.3750E-2, 3.7470E-2, 4.1528E-2, 4.6170E-2, 4.8708E-2, 5.1401E-2, 5.7297E-2, 6.3943E-2, 7.1441E-2}, 
  {0.04, 6.9619E-2, 7.6098E-2, 8.3249E-2, 9.1150E-2, 9.5406E-2, 9.9881E-2, 1.0954E-1, 1.2023E-1, 1.3208E-1}, 
  {0.05, 1.1522E-1, 1.2470E-1, 1.3504E-1, 1.4632E-1, 1.5234E-1, 1.5864E-1, 1.7211E-1, 1.8686E-1, 2.0304E-1}, 
  {0.06, 1.6931E-1, 1.8179E-1, 1.9530E-1, 2.0991E-1, 2.1767E-1, 2.2576E-1, 2.4295E-1, 2.6165E-1, 2.8201E-1}, 
  {0.08, 2.9540E-1, 3.1361E-1, 3.3312E-1, 3.5403E-1, 3.6506E-1, 3.7650E-1, 4.0067E-1, 4.2673E-1, 4.5488E-1}, 
  {0.1,  4.3613E-1, 4.5956E-1,  4.852E-1, 5.1115E-1, 5.2514E-1, 5.3961E-1, 5.7008E-1, 6.0277E-1, 6.3793E-1}, 
  {0.15, 8.1014E-1, 8.4453E-1, 8.8093E-1, 9.1954E-1, 9.3973E-1, 9.6056E-1, 1.0043, 1.0509, 1.1008}, 
  {0.2, 1.1888, 1.2319, 1.2774, 1.3255, 1.3506, 1.3765, 1.4308, 1.4886, 1.5504}, 
  {0.3, 1.8422, 1.8983, 1.9575, 2.0201, 2.0528, 2.0864, 2.1569, 2.2319, 2.3120}, 
  {0.4, 2.3984, 2.4642, 2.5336, 2.6070, 2.6452, 2.6847, 2.7674, 2.8554, 2.9494}, 
  {0.5, 2.8181, 2.8915, 2.9690, 3.0509, 3.0937, 3.1378, 3.2301, 3.3285, 3.4337}, 
  {0.6, 3.1698, 3.2494, 3.3336, 3.4226, 3.4691, 3.5171, 3.6175, 3.7246, 3.8391}, 
  {0.7, 3.4652, 3.5502, 3.6400, 3.7351, 3.7848, 3.8360, 3.9433, 4.0578, 4.1804}, 
  {0.8, 3.7392, 3.8289, 3.9236, 4.0239, 4.0764, 4.1304, 4.2438, 4.3648, 4.4944}, 
  {1.0, 4.2295, 4.3269, 4.4299, 4.5391, 4.5962, 4.6551, 4.7786, 4.9106, 5.0520}, 
  {1.2, 4.5777, 4.6814, 4.7912, 4.9076, 4.9685, 5.0314, 5.1633, 5.3043, 5.4555}, 
  {1.4, 4.9001, 5.0092, 5.1247, 5.2473, 5.3114, 5.3776, 5.5166, 5.6653, 5.8249}, 
  {1.5, 5.0434, 5.1550, 5.2731, 5.3984, 5.4640, 5.5317, 5.6739, 5.8260, 5.9893}, 
  {1.7, 5.3011, 5.4170, 5.5398, 5.6701, 5.7373, 5.8088, 5.9568, 6.1152, 6.2853}, 
  {2.0, 5.6308, 5.7523, 5.8811, 6.0174, 6.0896, 6.1636, 6.3192, 6.4857, 6.6647}, 
  {3.0, 6.4224, 6.5580, 6.7019, 6.8549, 6.9351, 7.0180, 7.1925, 7.3795, 7.5808}, 
  {5.0, 7.3338, 7.4872, 7.6500, 7.8235, 7.9146, 8.0087, 8.2071, 8.4200, 8.6496}, 
  {7.0, 7.7938, 7.9588, 8.1342, 8.3211, 8.4193, 8.5209, 8.7350, 8.9651, 9.2133}
  }; 
 
  static const G4double bls3[28][9] = { 
  {0.005, 1.2895E-3, 1.3670E-3, 1.4524E-3, 1.6524E-3, 1.9078E-3, 2.2414E-3, 2.6889E-3, 3.3006E-3}, 
  {0.007, 2.6467E-3, 2.8242E-3, 3.0238E-3, 3.5045E-3, 4.1260E-3, 4.9376E-3, 6.0050E-3, 7.4152E-3}, 
  {0.01, 6.1472E-3, 6.6086E-3, 7.1246E-3, 8.3491E-3, 9.8871E-3, 1.1822E-2, 1.4261E-2, 1.7335E-2}, 
  {0.015, 1.63316E-2, 1.7572E-2, 1.8932E-2, 2.2053E-2, 2.5803E-2, 3.0308E-2, 3.5728E-2, 4.2258E-2}, 
  {0.02, 3.2634E-2, 3.4900E-2, 3.7348E-2, 4.2850E-2, 4.9278E-2, 5.6798E-2, 6.5610E-2, 7.5963E-2}, 
  {0.03, 7.9907E-2, 8.4544E-2, 8.9486E-2, 1.0032E-1, 1.1260E-1, 1.2656E-1, 1.4248E-1, 1.6071E-1}, 
  {0.04, 1.4523E-1, 1.5237E-1, 1.5985E-1, 1.7614E-1, 1.9434E-1, 2.1473E-1, 2.3766E-1, 2.6357E-1}, 
  {0.05, 2.2082E-1, 2.3036E-1, 2.4038E-1, 2.6199E-1, 2.8590E-1, 3.1248E-1, 3.4212E-1, 3.7536E-1}, 
  {0.06, 3.0423E-1, 3.1611E-1, 3.2854E-1, 3.5522E-1, 3.8459E-1, 4.1704E-1, 4.5306E-1, 4.9326E-1}, 
  {0.08, 4.8536E-1, 5.0156E-1, 5.1846E-1, 5.5453E-1, 5.9397E-1, 6.3728E-1, 6.8507E-1, 7.3810E-1}, 
  {0.1, 6.7586E-1, 6.9596E-1, 7.1688E-1, 7.6141E-1, 8.0992E-1, 8.6301E-1, 9.2142E-1, 9.8604E-1}, 
  {0.15, 1.1544, 1.1828, 1.2122, 1.2746, 1.3424, 1.4163, 1.4974, 1.5868}, 
  {0.2, 1.6167, 1.6517, 1.6880, 1.7650, 1.8484, 1.9394, 2.0390, 2.1489}, 
  {0.3, 2.3979, 2.4432, 2.4902, 2.5899, 2.6980, 2.8159, 2.9451, 3.0876}, 
  {0.4, 3.0502, 3.1034, 3.1586, 3.2758, 3.4030, 3.5416, 3.6938, 3.8620}, 
  {0.5, 3.5466, 3.6062, 3.6681, 3.7994, 3.9421, 4.0978, 4.2688, 4.4580}, 
  {0.6, 3.9620, 4.0270, 4.0945, 4.2378, 4.3935, 4.5636, 4.7506, 4.9576}, 
  {0.7, 4.3020, 4.3715, 4.4438, 4.5974, 4.7644, 4.9470, 5.1478, 5.3703}, 
  {0.8, 4.6336, 4.7072, 4.7837, 4.9463, 5.1233, 5.3169, 5.5300, 5.7661}, 
  {1.0, 5.2041, 5.2845, 5.3682, 5.5462, 5.7400, 5.9523, 6.1863, 6.4458}, 
  {1.2, 5.6182, 5.7044, 5.7940, 5.9848, 6.1927, 6.4206, 6.6719, 6.9510}, 
  {1.4, 5.9967, 6.0876, 6.1823, 6.3839, 6.6038, 6.8451, 7.1113, 7.4071}, 
  {1.5, 6.1652, 6.2583, 6.3553, 6.5618, 6.7871, 7.0344, 7.3073, 7.6107}, 
  {1.7, 6.4686, 6.5657, 6.6668, 6.8823, 7.1175, 7.3757, 7.6609, 7.9782}, 
  {2.0, 6.8577, 6.9600, 7.0666, 7.2937, 7.5418, 7.8143, 8.1156, 8.4510}, 
  {3.0, 7.7981, 7.9134, 8.0336, 8.2901, 8.5708, 8.8796, 9.2214, 9.6027}, 
  {5.0, 8.8978, 9.0297, 9.1673, 9.4612, 9.7834, 10.1384, 10.5323, 10.9722}, 
  {7.0, 9.4819, 9.6248, 9.7739, 10.0926, 10.4423, 10.8282, 11.2565, 11.7356}
  }; 

  static const G4double bll1[28][10] = { 
  {0.005, 3.6324E-5, 4.0609E-5, 4.5430E-5, 5.6969E-5, 7.1625E-5, 9.0279E-5, 1.1407E-4, 1.2834E-4, 1.4447E-4}, 
  {0.007, 1.8110E-4, 2.0001E-4, 2.2099E-4, 2.7006E-4, 3.3049E-4, 4.0498E-4, 4.9688E-4, 5.5061E-4, 6.1032E-4}, 
  {0.01, 8.6524E-4, 9.4223E-4, 1.0262E-3, 1.2178E-3, 1.4459E-3, 1.7174E-3, 2.0405E-3, 2.2245E-3, 2.4252E-3}, 
  {0.015, 4.2293E-3, 4.5314E-3, 4.8551E-3, 5.5731E-3, 6.3968E-3, 7.3414E-3, 8.4242E-3, 9.0236E-3, 9.6652E-3}, 
  {0.02, 1.1485E-2, 1.2172E-2, 1.2900E-2, 1.4486E-2, 1.6264E-2, 1.8256E-2, 2.0487E-2, 2.1702E-2, 2.2989E-2}, 
  {0.03, 3.9471E-2, 4.1270E-2, 4.3149E-2, 4.7163E-2, 5.1543E-2, 5.6423E-2, 6.1540E-2, 6.4326E-2, 6.7237E-2}, 
  {0.04, 8.4454E-2, 8.7599E-2, 9.0860E-2, 9.7747E-2, 1.0516E-1, 1.1313E-1, 1.2171E-1, 1.2625E-1, 1.3096E-1}, 
  {0.05, 1.4339E-1, 1.4795E-1, 1.5266E-1, 1.6253E-1, 1.7306E-1, 1.8430E-1, 1.9630E-1, 2.0261E-1, 2.0924E-1}, 
  {0.06, 2.1304E-1, 2.1899E-1, 2.2512E-1, 2.3794E-1, 2.5153E-1, 2.6596E-1, 2.8130E-1, 2.8934E-1, 2.9763E-1}, 
  {0.08, 3.7382E-1, 3.8241E-1, 3.9122E-1, 4.0955E-1, 4.2888E-1, 4.4928E-1, 4.7086E-1, 4.8212E-1, 4.9371E-1}, 
  {0.1, 5.5056E-1, 5.6151E-1, 5.7273E-1, 5.9601E-1, 6.2049E-1, 6.4627E-1, 6.7346E-1, 6.8762E-1, 7.0218E-1}, 
  {0.15, 1.0066, 1.0224, 1.0386, 1.0721, 1.1073, 1.1443, 1.1832, 1.2035, 1.2243}, 
  {0.2, 1.4376, 1.4572, 1.4773, 1.5188, 1.5624, 1.6083, 1.6566, 1.6817, 1.7076}, 
  {0.3, 2.1712, 2.1964, 2.2223, 2.2758, 2.3322, 2.3915, 2.4542, 2.4869, 2.5205}, 
  {0.4, 2.7500, 2.7793, 2.8094, 2.8719, 2.9377, 3.0072, 3.0807, 3.1192, 3.1587}, 
  {0.5, 3.2033, 3.2359, 3.2693, 3.3389, 3.4122, 3.4898, 3.5721, 3.6151, 3.6595}, 
  {0.6, 3.6038, 3.6391, 3.6753, 3.7506, 3.8303, 3.9146, 4.0042, 4.0511, 4.0995}, 
  {0.7, 3.9106, 3.9482, 3.9867, 4.0670, 4.1520, 4.2421, 4.3380, 4.3882, 4.4401}, 
  {0.8, 4.1790, 4.2185, 4.2590, 4.3437, 4.4333, 4.5285, 4.6298, 4.6830, 4.7380}, 
  {1.0, 4.6344, 4.6772, 4.7212, 4.8131, 4.9106, 5.0144, 5.1250, 5.1831, 5.2432}, 
  {1.2, 4.9787, 5.0242, 5.0711, 5.1689, 5.2729, 5.3837, 5.5050, 5.5642, 5.6287}, 
  {1.4, 5.2688, 5.3166, 5.3658, 5.4687, 5.5782, 5.6950, 5.8198, 5.8855, 5.9554}, 
  {1.5, 5.3966, 5.4454, 5.4957, 5.6009, 5.7128, 5.8323, 5.9601, 6.0274, 6.0972}, 
  {1.7, 5.6255, 5.6762, 5.7284, 5.8377, 5.9541, 6.0785, 6.2116, 6.2818, 6.3546}, 
  {2.0, 5.9170, 5.9701, 6.0248, 6.1395, 6.2618, 6.3925, 6.5327, 6.6066, 6.6833}, 
  {3.0, 6.6210, 6.6801, 6.7411, 6.8692, 7.0062, 7.1529, 7.3107, 7.3941, 7.4807}, 
  {5.0, 7.4620, 7.5288, 7.5977, 7.7428, 7.8982, 8.0653, 8.2454, 8.3409, 8.4402}, 
  {7.0, 7.7362, 7.8079, 7.8821, 8.0383, 8.2061, 8.3866, 8.5816, 8.6850, 8.7927}
  }; 

  static const G4double bll2[28][10] = { 
  {0.005, 1.8339E-4, 2.3330E-4, 2.9738E-4, 3.7977E-4, 4.2945E-4, 4.8582E-4, 6.2244E-4, 7.9858E-4, 1.0258E-3}, 
  {0.007, 7.5042E-4, 9.2355E-4, 1.1375E-3, 1.4021E-3, 1.5570E-3, 1.7292E-3, 2.1335E-3, 2.6335E-3, 3.2515E-3}, 
  {0.01, 2.8829E-3, 3.4275E-3, 4.0758E-3, 4.8457E-3, 5.2839E-3, 5.7617E-3, 6.8504E-3, 8.1442E-3, 9.6816E-3}, 
  {0.015, 1.1087E-2, 1.2716E-2, 1.4581E-2, 1.6717E-2, 1.7898E-2, 1.9163E-2, 2.1964E-2, 2.5173E-2, 2.8851E-2}, 
  {0.02, 2.5786E-2, 2.8922E-2, 3.2435E-2, 3.6371E-2, 3.8514E-2, 4.0784E-2, 4.5733E-2, 5.1288E-2, 5.7531E-2}, 
  {0.03, 7.3461E-2, 8.0264E-2, 8.7705E-2, 9.5852E-2, 1.0021E-1, 1.0478E-1, 1.1458E-1, 1.2535E-1, 1.3721E-1}, 
  {0.04, 1.4094E-1, 1.5172E-1, 1.6336E-1, 1.7596E-1, 1.8265E-1, 1.8962E-1, 2.0445E-1, 2.2058E-1, 2.3818E-1}, 
  {0.05, 2.2289E-1, 2.3762E-1, 2.5344E-1, 2.7045E-1, 2.7944E-1, 2.8877E-1, 3.0855E-1, 3.2995E-1, 3.5318E-1}, 
  {0.06, 3.1503E-1, 3.3361E-1, 3.5346E-1, 3.7473E-1, 3.8594E-1, 3.9756E-1, 4.2212E-1, 4.4861E-1, 4.7727E-1}, 
  {0.08, 5.1793E-1, 5.4368E-1, 5.7109E-1, 6.0032E-1, 6.1569E-1, 6.3159E-1, 6.6512E-1, 7.0119E-1, 7.4012E-1}, 
  {0.1, 7.3258E-1, 7.6481E-1, 7.9907E-1, 8.3556E-1, 8.5472E-1, 8.7454E-1, 9.1630E-1, 9.6119E-1, 1.0096}, 
  {0.15, 1.2677, 1.3137, 1.3626, 1.4147, 1.4421, 1.4703, 1.5300, 1.5942, 1.6636}, 
  {0.2, 1.7615, 1.8188, 1.8797, 1.9446, 1.9788, 2.0142, 2.0889, 2.1694, 2.2567}, 
  {0.3, 2.5909, 2.6658, 2.7457, 2.8312, 2.8762, 2.9231, 3.0222, 3.1295, 3.2463}, 
  {0.4, 3.2417, 3.3302, 3.4249, 3.5266, 3.5803, 3.6361, 3.7546, 3.8835, 4.0242}, 
  {0.5, 3.7527, 3.8523, 3.9591, 4.0741, 4.1350, 4.1983, 4.3330, 4.4799, 4.6408}, 
  {0.6, 4.2013, 4.3103, 4.4274, 4.5537, 4.6206, 4.6904, 4.8390, 5.0013, 5.1796}, 
  {0.7, 4.5493, 4.6664, 4.7925, 4.9286, 5.0009, 5.0762, 5.2370, 5.4129, 5.6066}, 
  {0.8, 4.8537, 4.9780, 5.1119, 5.2568, 5.3338, 5.4141, 5.5857, 5.7738, 5.9811}, 
  {1.0, 5.3701, 5.5066, 5.6540, 5.8138, 5.8989, 5.9878, 6.1780, 6.3870, 6.6179}, 
  {1.2, 5.7648, 5.9114, 6.0701, 6.2424, 6.3343, 6.4303, 6.6361, 6.8626, 7.1137}, 
  {1.4, 6.0976, 6.2530, 6.4214, 6.6044, 6.7021, 6.8043, 7.0237, 7.2655, 7.5338}, 
  {1.5, 6.2447, 6.4041, 6.5768, 6.7647, 6.8650, 6.9700, 7.1954, 7.4442, 7.7203}, 
  {1.7, 6.5087, 6.6752, 6.8558, 7.0526, 7.1578, 7.2679, 7.5045, 7.7660, 8.0565}, 
  {2.0, 6.8458, 7.0218, 7.2129, 7.4213, 7.5328, 7.6496, 7.9010, 8.1791, 8.4886}, 
  {3.0, 7.6647, 7.8644, 8.0819, 8.3189, 8.4475, 8.5814, 8.8702, 9.1908, 9.5488}, 
  {5.0, 8.6515, 8.8816, 9.1330, 9.4090, 9.5572, 9.7132, 10.0504, 10.4259, 10.8466}, 
  {7.0, 9.0221, 9.2724, 9.5464, 9.8477, 10.0099, 10.1805, 10.5499, 10.9622, 11.4250}
  }; 

  static const G4double bll3[28][9] = { 
  {0.005, 1.3190E-3, 1.4961E-3, 1.6974E-3, 2.1858E-3, 2.8163E-3, 3.6302E-3, 4.6814E-3, 6.0395E-3}, 
  {0.007, 4.0158E-3, 4.4623E-3, 4.9592E-3, 6.1257E-3, 7.5675E-3, 9.3502E-3, 1.1556E-2, 1.4290E-2}, 
  {0.01, 1.1509E-2, 1.2548E-2, 1.3681E-2, 1.6263E-2, 1.9336E-2, 2.2999E-2, 2.7370E-2, 3.2603E-2}, 
  {0.015, 3.3070E-2, 3.5408E-2, 3.7914E-2, 4.3483E-2, 4.9898E-2, 5.7304E-2, 6.5884E-2, 7.5861E-2}, 
  {0.02, 6.4555E-2, 6.8394E-2, 7.2472E-2, 8.1413E-2, 9.1539E-2, 1.0304E-1, 1.1617E-1, 1.3121E-1}, 
  {0.03, 1.5030E-1, 1.5101E-1, 1.5844E-1, 1.7451E-1, 1.9244E-1, 2.1244E-1, 2.3496E-1, 2.6044E-1}, 
  {0.04, 2.5743E-1, 2.6774E-1, 2.7855E-1, 3.0180E-1, 3.2751E-1, 3.5608E-1, 3.8803E-1, 4.2401E-1}, 
  {0.05, 3.7846E-1, 3.9195E-1, 4.0607E-1, 4.3635E-1, 4.6973E-1, 5.0672E-1, 5.4798E-1, 5.9436E-1}, 
  {0.06, 5.0839E-1, 5.2497E-1, 5.4230E-1, 5.7943E-1, 6.2028E-1, 6.6549E-1, 7.1589E-1, 7.7252E-1}, 
  {0.08, 7.8230E-1, 8.0474E-1, 8.2818E-1, 8.7836E-1, 9.3355E-1, 9.9462E-1, 1.0627, 1.1394}, 
  {0.1, 1.0621, 1.0900, 1.1192, 1.1816, 1.2503, 1.3265, 1.4116, 1.5076}, 
  {0.15, 1.7389, 1.7790, 1.8210, 1.9112, 2.0108, 2.1217, 2.2462, 2.3876}, 
  {0.2, 2.3516, 2.4024, 2.4556, 2.5701, 2.6971, 2.8391, 2.9994, 3.1822}, 
  {0.3, 3.3741, 3.4427, 3.5148, 3.6706, 3.8445, 4.0404, 4.2631, 4.5193}, 
  {0.4, 4.1788, 4.2620, 4.3496, 4.5398, 4.7530, 4.9944, 5.2703, 5.5895}, 
  {0.5, 4.8180, 4.9137, 5.0146, 5.2341, 5.4811, 5.7618, 6.0840, 6.4583}, 
  {0.6, 5.3765, 5.4830, 5.5954, 5.8406, 6.1173, 6.4326, 6.7958, 7.2191}, 
  {0.7, 5.8208, 5.9369, 6.0596, 6.3275, 6.6306, 6.9769, 7.3767, 7.8440}, 
  {0.8, 6.2109, 6.3355, 6.4674, 6.7558, 7.0827, 7.4570, 7.8900, 8.3972}, 
  {1.0, 6.8747, 7.0142, 7.1621, 7.4861, 7.8546, 8.2778, 8.7690, 9.3464}, 
  {1.2, 7.3933, 7.5454, 7.7068, 8.0612, 8.4652, 8.9302, 9.4713, 10.1090}, 
  {1.4, 7.8331, 7.9967, 8.1694, 8.5502, 8.9851, 9.4866, 10.0713, 10.7619}, 
  {1.5, 8.0286, 8.1967, 8.3753, 8.7681, 9.2181, 9.7352, 10.3399, 11.0546}, 
  {1.7, 8.3813, 8.5585, 8.7469, 9.1618, 9.6367, 10.1856, 10.8270, 11.5863}, 
  {2.0, 8.8352, 9.0245, 9.2260, 9.6701, 10.1793, 10.7688, 11.4590, 12.2775}, 
  {3.0, 9.9511, 10.1714, 10.4062, 10.9254, 11.5229, 12.2172, 13.0332, 14.0048}, 
  {5.0, 11.3211, 11.5818, 11.8601, 12.4771, 13.1898, 14.0213, 15.0024, 16.1752}, 
  {7.0, 11.9480, 12.2357, 12.5432, 13.2260, 14.0164, 14.9404, 16.0330, 17.3420}
  };
                             
  G4double b, bs; 
  for(i=0; i<nEtaK; ++i) {

    G4double et = Eta[i];
    G4double loget = G4Log(et);

    for(j=0; j<nK; ++j) {

      if(j < 10) { b = bk2[i][10-j]; }
      else       { b = bk1[i][20-j]; }

      CK[j][i] = SK[j]*loget + TK[j] - b;

      if(i == nEtaK-1) { 
        ZK[j] = et*(et*et*CK[j][i] - et*UK[j] - VK[j]); 
        //G4cout << "i= " << i << " j= " << j 
        //       << " CK[j][i]= " <<  CK[j][i]
        //       << " ZK[j]= " << ZK[j] << "  b= " << b << G4endl;  
      } 
    }
    //G4cout << G4endl;
    if(i < nEtaL) {
      //G4cout << "  LShell:" <<G4endl;
      for(j=0; j<nL; ++j) {

        if(j < 8) {
          bs = bls3[i][8-j];
          b  = bll3[i][8-j];
        } else if(j < 17) {
          bs = bls2[i][17-j];
          b  = bll2[i][17-j];
        } else {
          bs = bls1[i][26-j];
          b  = bll1[i][26-j];
        }
        G4double c = SL[j]*loget + TL[j]; 
        CL[j][i] = c - bs - 3.0*b;
        if(i == nEtaL-1) { 
          VL[j] = et*(et*CL[j][i] - UL[j]); 
          //G4cout << "i= " << i << " j= " << j 
          //         << " CL[j][i]= " <<  CL[j][i]
          //         << " VL[j]= " << VL[j] << " b= " << b << " bs= " << bs 
          //         << " et= " << et << G4endl; 
          //" UL= " << UL[j] << " TL= " << TL[j] << " SL= " << SL[j] <<G4endl;  
        }
      }
      //G4cout << G4endl;
    }
  }

  static const G4double xzk[34] = { 11.7711,
    13.3669, 15.5762, 17.1715, 18.7667, 20.8523, 23.0606, 24.901, 26.9861, 29.4394, 31.77,
    34.3457, 37.4119, 40.3555, 42.3177, 44.7705, 47.2234, 50.78, 53.8458, 56.4214, 58.3834,
    60.9586, 63.6567, 66.5998, 68.807, 71.8728, 74.5706, 77.3911, 81.8056, 85.7297, 89.8988,
                             93.4549, 96.2753, 99.709};
  static const G4double yzk[34] = { 0.70663,
    0.72033, 0.73651, 0.74647, 0.75518, 0.76388, 0.77258, 0.78129, 0.78625, 0.7937, 0.79991,
    0.80611, 0.8123, 0.8185, 0.82097, 0.82467, 0.82838, 0.83457, 0.83702, 0.84198, 0.8432,
    0.84565, 0.84936, 0.85181, 0.85303, 0.85548, 0.85794, 0.8604, 0.86283, 0.86527, 0.86646,
                             0.86891, 0.87011, 0.87381};

  static const G4double xzl[36] = { 15.5102,
    16.7347, 17.9592, 19.551, 21.0204, 22.6122, 24.9388, 27.3878, 29.5918, 31.3061, 32.898,
    34.4898, 36.2041, 38.4082, 40.3674, 42.5714, 44.898, 47.4694, 49.9184, 52.7347, 55.9184,
    59.3469, 61.9184, 64.6122, 67.4286, 71.4694, 75.2653, 78.3265, 81.2653, 85.551, 88.7347,
                             91.551, 94.2449, 96.449, 98.4082, 99.7551};
  static const G4double yzl[36] = { 0.29875,
    0.31746, 0.33368, 0.35239, 0.36985, 0.38732, 0.41102, 0.43472, 0.45343, 0.4659, 0.47713,
    0.4896, 0.50083, 0.51331, 0.52328, 0.53077, 0.54075, 0.54823, 0.55572, 0.56445, 0.57193,
    0.58191, 0.5869, 0.59189, 0.60062, 0.60686, 0.61435, 0.61809, 0.62183, 0.62931, 0.6343,
                              0.6368, 0.64054, 0.64304, 0.64428, 0.64678};

  ThetaK = new G4LPhysicsFreeVector(34, xzk[0], xzk[33]);
  ThetaL = new G4LPhysicsFreeVector(36, xzl[0], xzl[35]);
  for(i=0; i<34; ++i) { ThetaK->PutValues(i, xzk[i], yzk[i]); }
  for(i=0; i<36; ++i) { ThetaL->PutValues(i, xzl[i], yzl[i]); }
  ThetaK->SetSpline(true);
  ThetaL->SetSpline(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

