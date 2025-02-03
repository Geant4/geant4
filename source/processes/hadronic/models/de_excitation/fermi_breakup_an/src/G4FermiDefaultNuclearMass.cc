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
// G4FermiBreakUpAN alternative de-excitation model
// by A. Novikov (January 2025)
//

#include "G4FermiDefaultNuclearMass.hh"

#include <CLHEP/Units/PhysicalConstants.h>

using CLHEP::MeV;

G4FermiDefaultNuclearMass::G4FermiDefaultNuclearMass()
{
#define FERMI_ADD_NUCLEI(A, Z, MASS) emplace_back(G4FermiNucleiData{A, Z}, MASS)

  FERMI_ADD_NUCLEI(1_m, 0_c, 939.565 * MeV);
  FERMI_ADD_NUCLEI(1_m, 1_c, 938.272 * MeV);
  FERMI_ADD_NUCLEI(2_m, 1_c, 1875.61 * MeV);
  FERMI_ADD_NUCLEI(3_m, 1_c, 2808.92 * MeV);
  FERMI_ADD_NUCLEI(3_m, 2_c, 2808.39 * MeV);
  FERMI_ADD_NUCLEI(3_m, 3_c, 2821.62 * MeV);
  FERMI_ADD_NUCLEI(4_m, 1_c, 3750.09 * MeV);
  FERMI_ADD_NUCLEI(4_m, 2_c, 3727.38 * MeV);
  FERMI_ADD_NUCLEI(4_m, 3_c, 3749.77 * MeV);
  FERMI_ADD_NUCLEI(5_m, 1_c, 4689.85 * MeV);
  FERMI_ADD_NUCLEI(5_m, 2_c, 4667.68 * MeV);
  FERMI_ADD_NUCLEI(5_m, 3_c, 4667.62 * MeV);
  FERMI_ADD_NUCLEI(5_m, 4_c, 4692.57 * MeV);
  FERMI_ADD_NUCLEI(6_m, 1_c, 5630.33 * MeV);
  FERMI_ADD_NUCLEI(6_m, 2_c, 5605.53 * MeV);
  FERMI_ADD_NUCLEI(6_m, 3_c, 5601.52 * MeV);
  FERMI_ADD_NUCLEI(6_m, 4_c, 5605.3 * MeV);
  FERMI_ADD_NUCLEI(6_m, 5_c, 5633.73 * MeV);
  FERMI_ADD_NUCLEI(7_m, 1_c, 6569.08 * MeV);
  FERMI_ADD_NUCLEI(7_m, 2_c, 6545.51 * MeV);
  FERMI_ADD_NUCLEI(7_m, 3_c, 6533.83 * MeV);
  FERMI_ADD_NUCLEI(7_m, 4_c, 6534.18 * MeV);
  FERMI_ADD_NUCLEI(7_m, 5_c, 6545.58 * MeV);
  FERMI_ADD_NUCLEI(8_m, 2_c, 7482.54 * MeV);
  FERMI_ADD_NUCLEI(8_m, 3_c, 7471.37 * MeV);
  FERMI_ADD_NUCLEI(8_m, 4_c, 7454.85 * MeV);
  FERMI_ADD_NUCLEI(8_m, 5_c, 7472.32 * MeV);
  FERMI_ADD_NUCLEI(8_m, 6_c, 7483.95 * MeV);
  FERMI_ADD_NUCLEI(9_m, 2_c, 8423.36 * MeV);
  FERMI_ADD_NUCLEI(9_m, 3_c, 8406.87 * MeV);
  FERMI_ADD_NUCLEI(9_m, 4_c, 8392.75 * MeV);
  FERMI_ADD_NUCLEI(9_m, 5_c, 8393.31 * MeV);
  FERMI_ADD_NUCLEI(9_m, 6_c, 8409.29 * MeV);
  FERMI_ADD_NUCLEI(10_m, 2_c, 9363.09 * MeV);
  FERMI_ADD_NUCLEI(10_m, 3_c, 9346.46 * MeV);
  FERMI_ADD_NUCLEI(10_m, 4_c, 9325.5 * MeV);
  FERMI_ADD_NUCLEI(10_m, 5_c, 9324.44 * MeV);
  FERMI_ADD_NUCLEI(10_m, 6_c, 9327.57 * MeV);
  FERMI_ADD_NUCLEI(10_m, 7_c, 9350.16 * MeV);
  FERMI_ADD_NUCLEI(11_m, 3_c, 10285.6 * MeV);
  FERMI_ADD_NUCLEI(11_m, 4_c, 10264.6 * MeV);
  FERMI_ADD_NUCLEI(11_m, 5_c, 10252.5 * MeV);
  FERMI_ADD_NUCLEI(11_m, 6_c, 10254.0 * MeV);
  FERMI_ADD_NUCLEI(11_m, 7_c, 10267.2 * MeV);
  FERMI_ADD_NUCLEI(12_m, 3_c, 11225.3 * MeV);
  FERMI_ADD_NUCLEI(12_m, 4_c, 11201.0 * MeV);
  FERMI_ADD_NUCLEI(12_m, 5_c, 11188.7 * MeV);
  FERMI_ADD_NUCLEI(12_m, 6_c, 11174.9 * MeV);
  FERMI_ADD_NUCLEI(12_m, 7_c, 11191.7 * MeV);
  FERMI_ADD_NUCLEI(12_m, 8_c, 11205.8 * MeV);
  FERMI_ADD_NUCLEI(13_m, 3_c, 12166.2 * MeV);
  FERMI_ADD_NUCLEI(13_m, 4_c, 12141.0 * MeV);
  FERMI_ADD_NUCLEI(13_m, 5_c, 12123.4 * MeV);
  FERMI_ADD_NUCLEI(13_m, 6_c, 12109.5 * MeV);
  FERMI_ADD_NUCLEI(13_m, 7_c, 12111.2 * MeV);
  FERMI_ADD_NUCLEI(13_m, 8_c, 12128.5 * MeV);
  FERMI_ADD_NUCLEI(14_m, 4_c, 13078.8 * MeV);
  FERMI_ADD_NUCLEI(14_m, 5_c, 13062.0 * MeV);
  FERMI_ADD_NUCLEI(14_m, 6_c, 13040.9 * MeV);
  FERMI_ADD_NUCLEI(14_m, 7_c, 13040.2 * MeV);
  FERMI_ADD_NUCLEI(14_m, 8_c, 13044.8 * MeV);
  FERMI_ADD_NUCLEI(14_m, 9_c, 13068.3 * MeV);
  FERMI_ADD_NUCLEI(15_m, 4_c, 14020.1 * MeV);
  FERMI_ADD_NUCLEI(15_m, 5_c, 13998.8 * MeV);
  FERMI_ADD_NUCLEI(15_m, 6_c, 13979.2 * MeV);
  FERMI_ADD_NUCLEI(15_m, 7_c, 13968.9 * MeV);
  FERMI_ADD_NUCLEI(15_m, 8_c, 13971.2 * MeV);
  FERMI_ADD_NUCLEI(15_m, 9_c, 13984.6 * MeV);
  FERMI_ADD_NUCLEI(16_m, 4_c, 14959.3 * MeV);
  FERMI_ADD_NUCLEI(16_m, 5_c, 14938.5 * MeV);
  FERMI_ADD_NUCLEI(16_m, 6_c, 14914.5 * MeV);
  FERMI_ADD_NUCLEI(16_m, 7_c, 14906.0 * MeV);
  FERMI_ADD_NUCLEI(16_m, 8_c, 14895.1 * MeV);
  FERMI_ADD_NUCLEI(16_m, 9_c, 14910.0 * MeV);
  FERMI_ADD_NUCLEI(16_m, 10_c, 14922.8 * MeV);
  FERMI_ADD_NUCLEI(17_m, 5_c, 15876.6 * MeV);
  FERMI_ADD_NUCLEI(17_m, 6_c, 15853.4 * MeV);
  FERMI_ADD_NUCLEI(17_m, 7_c, 15839.7 * MeV);
  FERMI_ADD_NUCLEI(17_m, 8_c, 15830.5 * MeV);
  FERMI_ADD_NUCLEI(17_m, 9_c, 15832.8 * MeV);
  FERMI_ADD_NUCLEI(17_m, 10_c, 15846.8 * MeV);
  FERMI_ADD_NUCLEI(18_m, 5_c, 16816.2 * MeV);
  FERMI_ADD_NUCLEI(18_m, 6_c, 16788.7 * MeV);
  FERMI_ADD_NUCLEI(18_m, 7_c, 16776.4 * MeV);
  FERMI_ADD_NUCLEI(18_m, 8_c, 16762.0 * MeV);
  FERMI_ADD_NUCLEI(18_m, 9_c, 16763.2 * MeV);
  FERMI_ADD_NUCLEI(18_m, 10_c, 16767.1 * MeV);
  FERMI_ADD_NUCLEI(18_m, 11_c, 16786.3 * MeV);
  FERMI_ADD_NUCLEI(19_m, 5_c, 17754.6 * MeV);
  FERMI_ADD_NUCLEI(19_m, 6_c, 17727.7 * MeV);
  FERMI_ADD_NUCLEI(19_m, 7_c, 17710.7 * MeV);
  FERMI_ADD_NUCLEI(19_m, 8_c, 17697.6 * MeV);
  FERMI_ADD_NUCLEI(19_m, 9_c, 17692.3 * MeV);
  FERMI_ADD_NUCLEI(19_m, 10_c, 17695.0 * MeV);
  FERMI_ADD_NUCLEI(19_m, 11_c, 17705.7 * MeV);
  FERMI_ADD_NUCLEI(19_m, 12_c, 17724.1 * MeV);
  FERMI_ADD_NUCLEI(20_m, 5_c, 18694.5 * MeV);
  FERMI_ADD_NUCLEI(20_m, 6_c, 18664.4 * MeV);
  FERMI_ADD_NUCLEI(20_m, 7_c, 18648.1 * MeV);
  FERMI_ADD_NUCLEI(20_m, 8_c, 18629.6 * MeV);
  FERMI_ADD_NUCLEI(20_m, 9_c, 18625.3 * MeV);
  FERMI_ADD_NUCLEI(20_m, 10_c, 18617.7 * MeV);
  FERMI_ADD_NUCLEI(20_m, 11_c, 18631.1 * MeV);
  FERMI_ADD_NUCLEI(20_m, 12_c, 18641.3 * MeV);

#undef FERMI_ADD_NUCLEI
}
