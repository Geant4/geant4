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
// Class Description
// Name list of Elements for a high precision (based on evaluated data
// libraries) description of themal neutron scattering below 4 eV;
// Based on Thermal neutron scattering files
// from the evaluated nuclear data files ENDF/B-VI, Release2
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.
// Class Description - End

// 15-Nov-06 First implementation is done by T. Koi (SLAC/SCCS)
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//

#include "G4ParticleHPThermalScatteringNames.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
//#include "G4ParticleHPData.hh"

G4ParticleHPThermalScatteringNames::G4ParticleHPThermalScatteringNames()
{
	// --------------------------------------------------------------------------------------------------------------------------
	// Old Geant4 naming - before 23/03/2022 - TSL linked to ENDF/BVII.1 nuclear cross-section G4NDL4.5
	// --------------------------------------------------------------------------------------------------------------------------
   /*names.insert ( std::pair < G4String , G4String > ( "TS_Aluminium_Metal" , "al_metal" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_Beryllium_Metal" , "be_metal" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_Be_of_Beryllium_Oxide" , "be_beo" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_C_of_Graphite" , "graphite" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_D_of_Heavy_Water" , "d_heavy_water" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Water" , "h_water" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Zirconium_Hydride" , "h_zrh" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Polyethylene" , "h_polyethylene" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_Iron_Metal" , "fe_metal" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_O_of_Uranium_Dioxide" , "o_uo2" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_O_of_Beryllium_Oxide" , "o_beo" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_U_of_Uranium_Dioxide" , "u_uo2" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_U235_of_Uranium_Dioxide" , "u235_uo2" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_U238_of_Uranium_Dioxide" , "u238_uo2" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_Zr_of_Zirconium_Hydride" , "zr_zrh" ) ); //// ENDF-B71
   

   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Para_Hydrogen" , "h_para_h2" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Ortho_Hydrogen" , "h_ortho_h2" ) ); 
   
   names.insert ( std::pair < G4String , G4String > ( "TS_D_of_Para_Deuterium" , "d_para_d2" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_D_of_Ortho_Deuterium" , "d_ortho_d2" ) ); 
   
   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Liquid_Methane", "h_l_ch4" ) ); 
   names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Solid_Methane", "h_s_ch4" ) ); */
   

	// --------------------------------------------------------------------------------------------------------------------------
	// New Geant4 naming - after 23/03/2022 - TSL linked to JEFF-3.3 nuclear cross-section G4NDL4.6
	// --------------------------------------------------------------------------------------------------------------------------
	///23/03/2022 - Added by L. Thulliez (CEA-Saclay)
	names.insert ( std::pair < G4String , G4String > ( "TS_Benzene", "h_benzen" ) ); 		 ///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Para_Hydrogen", "h_para_h2" ) );  ///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_D_of_Para_Deuterium", "d_para_d2" ) ); ///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Ortho_Hydrogen", "h_ortho_h2" ) ); ///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_D_of_Ortho_Deuterium", "d_ortho_d2" ) );///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_O_of_Uranium_Dioxide", "o_uo2" ) ); 		 ///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_O_of_Ice", "o_ice" ) ); 		             ///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_O_of_Heavy_Water", "o_heavy_water" ) ); 	 ///ENDF/BVIII.0 and JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_O_of_Beryllium_Oxide", "o_beo" ) ); 	     ///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_N_of_Uranium_Nitride", "n_un" ) ); 	     ///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Liquid_Methane", "h_l_ch4" ) );      ///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Zirconium_Hydride", "h_zrh" ) ); ///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Yttrium_Hydride", "h_yh2" ) ); 		 ///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Ice", "h_ice" ) ); 		 ///ENDF/BVIII.0 and JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Water", "h_water" ) ); 	 ///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Polyethylene", "h_polyethylene" ) ); 	///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_PolymethylMethacrylate", "h_c5o2h8" ) ); ///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_D_of_Heavy_Water", "d_heavy_water" ) ); 		 ///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_C_of_Graphite", "graphite" ) ); ///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_C_of_Silicium_Carbide", "c_sic" ) ); 		            ///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_C_of_Graphite_Porosity_30percent", "graphite_30p" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_C_of_Graphite_Porosity_10percent", "graphite_10p" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_Beryllium_Metal", "be_metal")); ///ENDF/BVIII.0 and JEFFF.3.3 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Be_of_Beryllium_Oxide", "be_beo" ) ); 		///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Iron_Metal", "fe_metal" ) ); 		        ///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Zr90_of_Zirconium_Hydride", "zr90_zrh" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Zr91_of_Zirconium_Hydride", "zr91_zrh" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Zr92_of_Zirconium_Hydride", "zr92_zrh" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Zr94_of_Zirconium_Hydride", "zr94_zrh" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Zr96_of_Zirconium_Hydride", "zr96_zrh" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Y_of_Yttrium_Hydride", "y_yh2" ) ); 		    ///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_U235_of_Uranium_Dioxide", "u235_uo2" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_U238_of_Uranium_Dioxide", "u238_uo2" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_U235_of_Uranium_Nitride", "u235_un" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_U238_of_Uranium_Nitride", "u238_un" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si28_of_SiO2_beta", "si28_sio2_beta" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si29_of_SiO2_beta", "si29_sio2_beta" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si30_of_SiO2_beta", "si30_sio2_beta" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si28_of_SiO2_alpha", "si28_sio2_alpha" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si29_of_SiO2_alpha", "si29_sio2_alpha" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si30_of_SiO2_alpha", "si30_sio2_alpha" ) ); 	///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si28_of_Silicium_Carbide", "si28_sic" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si29_of_Silicium_Carbide", "si29_sic" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si30_of_Silicium_Carbide", "si30_sic" ) ); 	///ENDF/BVIII.0 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Solid_Methane", "h_s_ch4" ) ); 		 ///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Aluminium_Metal", "al_metal" ) ); 		 ///ENDF/BVIII.0 and ENDF/BVII.1 
	names.insert ( std::pair < G4String , G4String > ( "TS_Al_of_Sapphir_SingleCrystal", "al_al2o3_singlecrystal" ) ); 	///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_Ca_of_CaH2", "ca_cah2" ) ); 		 ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_CaH2", "h_cah2" ) ); 		 ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Mesitylene_phaseII", "h_mesitylene_phaseII" ) ); 	 ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_O_of_Sapphir_SingleCrystal", "o_al2o3_singlecrystal" ) ); ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_H_of_Toluene", "h_toluene" ) ); 		 ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si30_of_SiO2_SingleCrystal", "si30_sio2_singlecrystal" ) ); ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si29_of_SiO2_SingleCrystal", "si29_sio2_singlecrystal" ) ); ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_Si28_of_SiO2_SingleCrystal", "si28_sio2_singlecrystal" ) ); ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_Mg26_of_Magnesium_Metal", "mg26_magnesium" ) ); 		 ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_Mg25_of_Magnesium_Metal", "mg25_magnesium" ) ); 		 ///JEFFF.3.3 
	names.insert ( std::pair < G4String , G4String > ( "TS_Mg24_of_Magnesium_Metal", "mg24_magnesium" ) ); 		 ///JEFFF.3.3 

   nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_BERYLLIUM_OXIDE" , "Be" ) , "be_beo" ) );
   nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_BERYLLIUM_OXIDE" , "O" ) , "o_beo" ) );
   nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_GRAPHITE" , "C" ) , "graphite" ) );
   nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_POLYETHYLENE" , "H" ) , "h_polyethylene" ) );
   nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_URANIUM_OXIDE" , "O" ) , "o_uo2" ) );
   nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_URANIUM_OXIDE" , "U" ) , "u_uo2" ) );
   nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_WATER" , "H" ) , "h_water" ) );

   //nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_BENZENE" , "H" ) , "benzen" ) );
   //nist_names.insert ( std::pair < std::pair < G4String , G4String > , G4String > ( std::pair < G4String , G4String > ( "G4_BENZENE" , "C" ) , "benzen" ) );


	

}

G4ParticleHPThermalScatteringNames::~G4ParticleHPThermalScatteringNames()
{
;
}

G4bool G4ParticleHPThermalScatteringNames::IsThisThermalElement( G4String aname)
{
   G4bool result = false;
   if ( names.find ( aname ) != names.end() ) result = true; 
   return result;
}

G4bool G4ParticleHPThermalScatteringNames::IsThisThermalElement( G4String material , G4String element )
{
   G4bool result = false;
   if ( nist_names.find ( std::pair<G4String,G4String>(material,element) ) != nist_names.end() ) result = true; 
   return result;
}

                                                               //Name of G4Element , Name of NDL file
void G4ParticleHPThermalScatteringNames::AddThermalElement ( G4String nameG4Element , G4String filename)
{  
   if ( names.find ( nameG4Element ) == names.end() ) names.insert( std::pair<G4String,G4String>( nameG4Element , filename ) ); 
}
