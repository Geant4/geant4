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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// G4UCNMaterialPropertiesTable
// ----------------------------
//
// Derives from G4MaterialPropertiesTable and adds a double pointer to the
// MicroRoughnessTable
//
// This file includes the initialization and calculation of the mr-lookup
// tables. It also provides the functions to read from these tables and to
// calculate the probability for certain angles.
//
// For a closer description see the header file
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <fstream>

#include "G4UCNMaterialPropertiesTable.hh"
#include "G4UCNMicroRoughnessHelper.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4UCNMaterialPropertiesTable::G4UCNMaterialPropertiesTable()
                             : G4MaterialPropertiesTable() 
{
  theMicroRoughnessTable = nullptr;
  maxMicroRoughnessTable = nullptr;
  theMicroRoughnessTransTable = nullptr;
  maxMicroRoughnessTransTable = nullptr;

  theta_i_min =  0.*degree;
  theta_i_max = 90.*degree;

  Emin =    0.e-9*eV;
  Emax = 1000.e-9*eV;

  no_theta_i = 90;
  noE = 100;

  theta_i_step = (theta_i_max-theta_i_min)/(no_theta_i-1);
  E_step = (Emax-Emin)/(noE-1);

  b =  1*nm;
  w = 30*nm;

  AngCut = 0.01*degree;
}

G4UCNMaterialPropertiesTable::~G4UCNMaterialPropertiesTable()
{
  delete theMicroRoughnessTable;
  delete maxMicroRoughnessTable;
  delete theMicroRoughnessTransTable;
  delete maxMicroRoughnessTransTable;
}

G4double* G4UCNMaterialPropertiesTable::GetMicroRoughnessTable ()
{
  return theMicroRoughnessTable;
}

G4double* G4UCNMaterialPropertiesTable::GetMicroRoughnessTransTable ()
{
  return theMicroRoughnessTransTable;
}

void G4UCNMaterialPropertiesTable::
               LoadMicroRoughnessTables(G4double* pMicroRoughnessTable,
                                        G4double* pmaxMicroRoughnessTable,
                                        G4double* pMicroRoughnessTransTable,
                                        G4double* pmaxMicroRoughnessTransTable)
{
  theMicroRoughnessTable      = pMicroRoughnessTable;
  maxMicroRoughnessTable      = pmaxMicroRoughnessTable;
  theMicroRoughnessTransTable = pMicroRoughnessTransTable;
  maxMicroRoughnessTransTable = pmaxMicroRoughnessTransTable;
}

void G4UCNMaterialPropertiesTable::InitMicroRoughnessTables()
{
  G4int NEdim     = 0;
  G4int Nthetadim = 0;

  // Checks if the number of angles is available and stores it

  if(ConstPropertyExists("MR_NBTHETA"))
  {
    Nthetadim = G4int(GetConstProperty("MR_NBTHETA") + 0.1);
  }

  // Checks if the number of energies is available and stores it

  if(ConstPropertyExists("MR_NBE"))
  {
    NEdim = G4int(GetConstProperty("MR_NBE") + 0.1);
  }

  //G4cout << "thetadim: " << Nthetadim << " , Edim: " << NEdim << G4endl;

  // If both dimensions of the lookup-table are non-trivial:
  // delete old tables if existing and allocate memory for new tables

  if (Nthetadim*NEdim > 0) {
    delete theMicroRoughnessTable;
    theMicroRoughnessTable = new G4double[Nthetadim * NEdim];
    delete maxMicroRoughnessTable;
    maxMicroRoughnessTable = new G4double[Nthetadim * NEdim];
    delete theMicroRoughnessTransTable;
    theMicroRoughnessTransTable = new G4double[Nthetadim * NEdim];
    delete maxMicroRoughnessTransTable;
    maxMicroRoughnessTransTable = new G4double[Nthetadim * NEdim];
  }
}

void G4UCNMaterialPropertiesTable::ComputeMicroRoughnessTables()
{
// Reads the parameters for the mr-probability computation from the
// corresponding material properties and stores it in the appropriate
// variables

  b = GetConstProperty("MR_RRMS");
  G4double b2 = b*b;
  w = GetConstProperty("MR_CORRLEN");
  G4double w2 = w*w;

  no_theta_i = G4int(GetConstProperty("MR_NBTHETA")+0.1);
  //G4cout << "no_theta: " << no_theta_i << G4endl;
  noE = G4int(GetConstProperty("MR_NBE")+0.1);
  //G4cout << "noE:      " << noE << G4endl;

  theta_i_min = GetConstProperty("MR_THETAMIN");
  theta_i_max = GetConstProperty("MR_THETAMAX");
  Emin = GetConstProperty("MR_EMIN");
  Emax = GetConstProperty("MR_EMAX");
  G4int AngNoTheta = G4int(GetConstProperty("MR_ANGNOTHETA")+0.1);
  G4int AngNoPhi = G4int(GetConstProperty("MR_ANGNOPHI")+0.1);
  AngCut = GetConstProperty("MR_ANGCUT");

  // The Fermi potential was saved in neV by P.F.

  G4double fermipot = GetConstProperty("FERMIPOT")*(1.e-9*eV);

  //G4cout << "Fermipot: " << fermipot/(1.e-9*eV) << "neV" << G4endl;

  G4double theta_i, E;

  // Calculates the increment in theta_i in the lookup-table
  theta_i_step = (theta_i_max-theta_i_min)/(no_theta_i-1);

  //G4cout << "theta_i_step: " << theta_i_step << G4endl;

  // Calculates the increment in energy in the lookup-table
  E_step = (Emax-Emin)/(noE-1);

  // Runs the lookup-table memory allocation
  InitMicroRoughnessTables();

  G4int counter = 0;

  //G4cout << "Calculating Microroughnesstables...";

  // Writes the mr-lookup-tables to files for immediate control

  std::ofstream dateir("MRrefl.dat",std::ios::out);
  std::ofstream dateit("MRtrans.dat",std::ios::out);

  //G4cout << theMicroRoughnessTable << G4endl;

  for (theta_i=theta_i_min; theta_i<=theta_i_max+1e-6; theta_i+=theta_i_step) {
      // Calculation for each cell in the lookup-table
      for (E=Emin; E<=Emax; E+=E_step) {
          *(theMicroRoughnessTable+counter) =
                      G4UCNMicroRoughnessHelper::GetInstance() ->
                      IntIplus(E, fermipot, theta_i, AngNoTheta, AngNoPhi,
                               b2, w2, maxMicroRoughnessTable+counter, AngCut);

	  *(theMicroRoughnessTransTable+counter) =
                      G4UCNMicroRoughnessHelper::GetInstance() -> 
                      IntIminus(E, fermipot, theta_i, AngNoTheta, AngNoPhi,
                                b2, w2, maxMicroRoughnessTransTable+counter,
                                AngCut);

          dateir << *(theMicroRoughnessTable+counter)      << G4endl;
          dateit << *(theMicroRoughnessTransTable+counter) << G4endl;

          counter++;

          //G4cout << counter << G4endl;
      }
  }

  dateir.close();
  dateit.close();

  // Writes the mr-lookup-tables to files for immediate control

  std::ofstream dateic("MRcheck.dat",std::ios::out);
  std::ofstream dateimr("MRmaxrefl.dat",std::ios::out);
  std::ofstream dateimt("MRmaxtrans.dat",std::ios::out);

  for (theta_i=theta_i_min; theta_i<=theta_i_max+1e-6; theta_i+=theta_i_step) {
      for (E=Emin; E<=Emax; E+=E_step) {

          // tests the GetXXProbability functions by writing the entries
          // of the lookup tables to files

          dateic  << GetMRIntProbability(theta_i,E)      << G4endl;
          dateimr << GetMRMaxProbability(theta_i,E)      << G4endl;
          dateimt << GetMRMaxTransProbability(theta_i,E) << G4endl;
      }
  }

  dateic.close();
  dateimr.close();
  dateimt.close();
}

G4double G4UCNMaterialPropertiesTable::
                  GetMRIntProbability(G4double theta_i, G4double Energy)
{
  if(theMicroRoughnessTable == nullptr)
  {
    G4cout << "Do not have theMicroRoughnessTable" << G4endl;
    return 0.;
  }

  // if theta_i or energy outside the range for which the lookup table is
  // calculated, the probability is set to zero

  //G4cout << "theta_i: " << theta_i/degree << "degree"
  //       << " theta_i_min: " << theta_i_min/degree << "degree"
  //       << " theta_i_max: " << theta_i_max/degree << "degree"
  //       << " Energy: " << Energy/(1.e-9*eV) << "neV"
  //       << " Emin: " << Emin/(1.e-9*eV) << "neV"
  //       << " Emax: " << Emax/(1.e-9*eV) << "neV" << G4endl;

  if(theta_i < theta_i_min || theta_i > theta_i_max || Energy < Emin ||
     Energy > Emax)
  {
    return 0.;
  }

  // Determines the nearest cell in the lookup table which contains
  // the probability

  G4int theta_i_pos = G4int((theta_i-theta_i_min)/theta_i_step+0.5);
  G4int E_pos = G4int((Energy-Emin)/E_step+0.5);

  // lookup table is onedimensional (1 row), energy is in rows,
  // theta_i in columns

  //G4cout << "E_pos: " << E_pos << " theta_i_pos: " << theta_i_pos << G4endl;
  //G4cout << "Probability: " << *(theMicroRoughnessTable+E_pos+theta_i_pos*(noE-1)) << G4endl;

  return *(theMicroRoughnessTable+E_pos+theta_i_pos*(noE - 1));
}

G4double G4UCNMaterialPropertiesTable::
                  GetMRIntTransProbability(G4double theta_i, G4double Energy)
{
  if(theMicroRoughnessTransTable == nullptr)
  {
    return 0.;
  }

  // if theta_i or energy outside the range for which the lookup table
  // is calculated, the probability is set to zero

  if(theta_i < theta_i_min || theta_i > theta_i_max || Energy < Emin ||
     Energy > Emax)
  {
    return 0.;
  }

  // Determines the nearest cell in the lookup table which contains
  // the probability

  G4int theta_i_pos = G4int((theta_i-theta_i_min)/theta_i_step+0.5);
  G4int E_pos = G4int((Energy-Emin)/E_step+0.5);

  // lookup table is onedimensional (1 row), energy is in rows,
  // theta_i in columns

  return *(theMicroRoughnessTransTable+E_pos+theta_i_pos*(noE - 1));
}

G4double G4UCNMaterialPropertiesTable::
                  GetMRMaxProbability(G4double theta_i, G4double Energy)
{
  if(maxMicroRoughnessTable == nullptr)
  {
    return 0.;
  }

  // if theta_i or energy outside the range for which the lookup table
  // is calculated, the probability is set to zero

  if(theta_i < theta_i_min || theta_i > theta_i_max || Energy < Emin ||
     Energy > Emax)
  {
    return 0.;
  }

  // Determines the nearest cell in the lookup table which contains
  // the probability

  G4int theta_i_pos = G4int((theta_i-theta_i_min)/theta_i_step+0.5);
  G4int E_pos = G4int((Energy-Emin)/E_step+0.5);

  // lookup table is onedimensional (1 row), energy is in rows,
  // theta_i in columns

  return *(maxMicroRoughnessTable+E_pos+theta_i_pos*noE);
}

void G4UCNMaterialPropertiesTable::
       SetMRMaxProbability(G4double theta_i, G4double Energy, G4double value)
{
  if(maxMicroRoughnessTable != nullptr)
  {
    if(theta_i < theta_i_min || theta_i > theta_i_max || Energy < Emin ||
       Energy > Emax)
    {}
    else
    {
      // Determines the nearest cell in the lookup table which contains
      // the probability

      G4int theta_i_pos = G4int((theta_i - theta_i_min) / theta_i_step + 0.5);
      G4int E_pos       = G4int((Energy - Emin) / E_step + 0.5);

      // lookup table is onedimensional (1 row), energy is in rows,
      // theta_i in columns

      *(maxMicroRoughnessTable + E_pos + theta_i_pos * noE) = value;
    }
  }
}

G4double G4UCNMaterialPropertiesTable::
                  GetMRMaxTransProbability(G4double theta_i, G4double Energy)
{
  if(maxMicroRoughnessTransTable == nullptr)
  {
    return 0.;
  }

  // if theta_i or energy outside the range for which the lookup table
  // is calculated, the probability is set to zero

  if(theta_i < theta_i_min || theta_i > theta_i_max || Energy < Emin ||
     Energy > Emax)
  {
    return 0.;
  }

  // Determines the nearest cell in the lookup table which contains
  // the probability

  G4int theta_i_pos = G4int((theta_i-theta_i_min)/theta_i_step+0.5);
  G4int E_pos = G4int((Energy-Emin)/E_step+0.5);

  // lookup table is onedimensional (1 row), energy is in rows,
  // theta_i in columns

  return *(maxMicroRoughnessTransTable+E_pos+theta_i_pos*noE);
}

void G4UCNMaterialPropertiesTable::
     SetMRMaxTransProbability(G4double theta_i, G4double Energy, G4double value)
{
  if(maxMicroRoughnessTransTable != nullptr)
  {
    if(theta_i < theta_i_min || theta_i > theta_i_max || Energy < Emin ||
       Energy > Emax)
    {}
    else
    {
      // Determines the nearest cell in the lookup table which contains
      // the probability

      G4int theta_i_pos = G4int((theta_i - theta_i_min) / theta_i_step + 0.5);
      G4int E_pos       = G4int((Energy - Emin) / E_step + 0.5);

      // lookup table is onedimensional (1 row), energy is in rows,
      // theta_i in columns

      *(maxMicroRoughnessTransTable + E_pos + theta_i_pos * noE) = value;
    }
  }
}

G4double G4UCNMaterialPropertiesTable::
                  GetMRProbability(G4double theta_i, G4double Energy,
                                   G4double fermipot,
                                   G4double theta_o, G4double phi_o)
{
  return G4UCNMicroRoughnessHelper::GetInstance()->
          ProbIplus(Energy, fermipot, theta_i, theta_o, phi_o, b, w, AngCut);
}

G4double G4UCNMaterialPropertiesTable::
                  GetMRTransProbability(G4double theta_i, G4double Energy,
                                        G4double fermipot,
                                        G4double theta_o, G4double phi_o)
{
  return G4UCNMicroRoughnessHelper::GetInstance()->
          ProbIminus(Energy, fermipot,theta_i, theta_o, phi_o, b, w, AngCut);
}

G4bool G4UCNMaterialPropertiesTable::ConditionsValid(G4double E,
                                                     G4double VFermi,
                                                     G4double theta_i)
{
  G4double k =   std::sqrt(2*neutron_mass_c2*E      / hbarc_squared);
  G4double k_l = std::sqrt(2*neutron_mass_c2*VFermi / hbarc_squared);

  //G4cout << " Energy: " << E/(1.e-9*eV) << "neV"
  //       << " VFermi: " << VFermi/(1.e-9*eV) << "neV"
  //       << " theta:  " << theta_i/degree << "degree" << G4endl;

  //G4cout << " Conditions - 2*b*k*cos(theta): " << 2*b*k*cos(theta_i)
  //       << ", 2*b*kl: "                       << 2*b*k_l << G4endl;

  // see eq. 17 of the Steyerl paper

  return 2 * b * k * std::cos(theta_i) < 1 && 2 * b * k_l < 1;
}

G4bool G4UCNMaterialPropertiesTable::TransConditionsValid(G4double E,
                                                          G4double VFermi,
                                                          G4double theta_i)
{
  G4double k2   = 2*neutron_mass_c2*E      / hbarc_squared;
  G4double k_l2 = 2*neutron_mass_c2*VFermi / hbarc_squared;

  if(E * (std::cos(theta_i) * std::cos(theta_i)) < VFermi)
  {
    return false;
  }

  G4double kS2 = k_l2 - k2;

  //G4cout << "Conditions; 2bk' cos(theta): " << 2*b*sqrt(kS2)*cos(theta_i) << 
  //          ", 2bk_l: " << 2*b*sqrt(k_l2) << G4endl;

  // see eq. 18 of the Steyerl paper

  return 2 * b * std::sqrt(kS2) * std::cos(theta_i) < 1 &&
         2 * b * std::sqrt(k_l2) < 1;
}

void G4UCNMaterialPropertiesTable::
       SetMicroRoughnessParameters(G4double ww, G4double bb,
                                   G4int no_theta, G4int no_E,
                                   G4double theta_min, G4double theta_max,
                                   G4double E_min, G4double E_max,
                                   G4int AngNoTheta, G4int AngNoPhi,
                                   G4double AngularCut)
{
  //G4cout << "Setting Microroughness Parameters...";

  // Removes an existing RMS roughness
  if(ConstPropertyExists("MR_RRMS"))
  {
    RemoveConstProperty("MR_RRMS");
  }

  // Adds a new RMS roughness
  AddConstProperty("MR_RRMS", bb);

  //G4cout << "b: " << bb << G4endl;

  // Removes an existing correlation length
  if(ConstPropertyExists("MR_CORRLEN"))
  {
    RemoveConstProperty("MR_CORRLEN");
  }

  // Adds a new correlation length
  AddConstProperty("MR_CORRLEN", ww);

  //G4cout << "w: " << ww << G4endl;

  // Removes an existing number of thetas
  if(ConstPropertyExists("MR_NBTHETA"))
  {
    RemoveConstProperty("MR_NBTHETA");
  }

  // Adds a new number of thetas
  AddConstProperty("MR_NBTHETA", (G4double)no_theta);

  //G4cout << "no_theta: " << no_theta << G4endl;

  // Removes an existing number of Energies
  if(ConstPropertyExists("MR_NBE"))
  {
    RemoveConstProperty("MR_NBE");
  }

  // Adds a new number of Energies
  AddConstProperty("MR_NBE", (G4double)no_E);

  //G4cout << "no_E: " << no_E << G4endl;

  // Removes an existing minimum theta
  if(ConstPropertyExists("MR_THETAMIN"))
  {
    RemoveConstProperty("MR_THETAMIN");
  }

  // Adds a new minimum theta
  AddConstProperty("MR_THETAMIN", theta_min);

  //G4cout << "theta_min: " << theta_min << G4endl;

  // Removes an existing maximum theta
  if(ConstPropertyExists("MR_THETAMAX"))
  {
    RemoveConstProperty("MR_THETAMAX");
  }

  // Adds a new maximum theta
  AddConstProperty("MR_THETAMAX", theta_max);

  //G4cout << "theta_max: " << theta_max << G4endl;

  // Removes an existing minimum energy
  if(ConstPropertyExists("MR_EMIN"))
  {
    RemoveConstProperty("MR_EMIN");
  }

  // Adds a new minimum energy
  AddConstProperty("MR_EMIN", E_min);

  //G4cout << "Emin: " << E_min << G4endl;

  // Removes an existing maximum energy
  if(ConstPropertyExists("MR_EMAX"))
  {
    RemoveConstProperty("MR_EMAX");
  }

  // Adds a new maximum energy
  AddConstProperty("MR_EMAX", E_max);

  //G4cout << "Emax: " << E_max << G4endl;

  // Removes an existing Theta angle number
  if(ConstPropertyExists("MR_ANGNOTHETA"))
  {
    RemoveConstProperty("MR_ANGNOTHETA");
  }

  // Adds a new Theta angle number
  AddConstProperty("MR_ANGNOTHETA", (G4double)AngNoTheta);

  //G4cout << "AngNoTheta: " << AngNoTheta << G4endl;

  // Removes an existing Phi angle number
  if(ConstPropertyExists("MR_ANGNOPHI"))
  {
    RemoveConstProperty("MR_ANGNOPHI");
  }

  // Adds a new Phi angle number
  AddConstProperty("MR_ANGNOPHI", (G4double)AngNoPhi);

  //G4cout << "AngNoPhi: " << AngNoPhi << G4endl;

  // Removes an existing angular cut
  if(ConstPropertyExists("MR_ANGCUT"))
  {
    RemoveConstProperty("MR_ANGCUT");
  }

  // Adds a new angle number
  AddConstProperty("MR_ANGCUT", AngularCut);

  //G4cout << "AngularCut: " << AngularCut/degree << "degree" << G4endl;

  // Starts the lookup table calculation

  ComputeMicroRoughnessTables();

  //G4cout << " *** DONE! ***" << G4endl;
}
