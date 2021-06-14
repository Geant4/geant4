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
/// \file hadronic/Hadr02/src/HadronicInelasticModelCRMC.cc
/// \brief Implementation of the HadronicInelasticModelCRMC class methods
//
//
//---------------------------------------------------------------------------
//
#ifdef G4_USE_CRMC

#include "HadronicInelasticModelCRMC.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4NucleiProperties.hh"

#include "Randomize.hh"

#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>

#define MAX_ENERGY_LAB_GEV 10000000.
#define MAX_ENERGY_CMS_GEV 30000.      // assuming that the target is <=100 times heavier than the projectile

#define IGNORE_PARTICLE_UNKNOWN_PDGID    false
#define USE_ENERGY_CORR                  false
#define ENERGY_NON_CONSERVATION_RESAMPLE false
#define ENERGY_NON_CONSERVATION_EMAX_GEV 0.999
#define ENERGY_NON_CONSERVATION_FRACTION_MAX 0.00001
#define ENERGY_NON_CONSERVATION_FRACTION_MAX_ATTEMPT 10
#define ENERGY_NON_CONSERVATION_FRACTION_MAX_ENERGYTRY ENERGY_NON_CONSERVATION_EMAX_GEV

#define SPLIT_MULTI_NEUTRONS_MAXN 10 
#define PARTICLE_MULTI_NEUTRONS_ERRORCODE -1 

#define ERROR_REPORT_EMAIL "andrii.tykhonov@SPAMNOTcern.ch"
#define CRMC_CONFIG_FILE_ENV_VARIABLE "CRMC_CONFIG_FILE"

//***********************************
// CRMC ION DEFINITION
// ID = CRMC_ION_COEF_0 + 
//      CRMC_ION_COEF_Z * Z + 
//      CRMC_ION_COEF_A * A
//
#define CRMC_ION_COEF_0  1000000000
#define CRMC_ION_COEF_Z  10000
#define CRMC_ION_COEF_A  10
//***********************************

HadronicInelasticModelCRMC::HadronicInelasticModelCRMC(int model, const G4String& modelName):
  G4HadronicInteraction(modelName), fPrintDebug(false)
{
	SetMaxEnergy(MAX_ENERGY_LAB_GEV * GeV);



	//int model = 1; // Epos (use temporary), it is faster
	//int model = 12; // Dpmjet
	int seed = 123456789; 
	//int seed = CLHEP::HepRandom::getTheSeed();  // Returns 0 which is invalid
	int produce_tables     = 0;                         // CRMC default, see CRMCoptions.cc in the CRMC package
	fTypeOutput            = 0;                         // CRMC default, see CRMCoptions.cc in the CRMC package
	static std::string crmc_param = GetCrmcParamPath(); //"crmc.param"; // CRMC default, see CRMCoptions.cc in the CRMC package
	
	fInterface = new CRMCinterface();
	fInterface->init(model);

	// open FORTRAN IO at first call
	fInterface->crmc_init(MAX_ENERGY_CMS_GEV,seed,model,produce_tables,fTypeOutput,crmc_param.c_str(),"",0);

	// final state
	finalState = new G4HadFinalState();

	// geant4 particle helpers (tables)
	fParticleTable = G4ParticleTable::GetParticleTable();
	fIonTable      = fParticleTable->GetIonTable();

}

std::string HadronicInelasticModelCRMC::GetCrmcParamPath(){
	std::string crmcParamPath = std::getenv(CRMC_CONFIG_FILE_ENV_VARIABLE);
	if (crmcParamPath==""){
		std::ostringstream errorstr;
		errorstr<<"CRMC ERROR: could not find crmc param file, please check "<< CRMC_CONFIG_FILE_ENV_VARIABLE <<" envornoment variable!";
		std::string error(errorstr.str());
		std::cout<<error<<std::endl;
		throw error;
	}
	std::cout<< "Using CRMC parameter file: " << crmcParamPath << std::endl;
	return crmcParamPath;
}

HadronicInelasticModelCRMC::~HadronicInelasticModelCRMC()
{}

G4HadFinalState * HadronicInelasticModelCRMC::ApplyYourself(const G4HadProjectile &aTrack, G4Nucleus &targetNucleus){

	//* leanup data vectors
	gCRMC_data.Clean();

	//* cleanup geant4 final state vector
	finalState->Clear();
	finalState->SetStatusChange(G4HadFinalStateStatus::stopAndKill); // TODO: check: inelastic collisions kills previos particles?

	//* git input particles parameters
	int id_proj  = aTrack.GetDefinition()->GetPDGEncoding();
	int id_targ  = targetNucleus.GetZ_asInt()*10000 +  targetNucleus.GetA_asInt()*10;
	double p_proj = aTrack.Get4Momentum().pz() / GeV; 
	double e_proj = aTrack.Get4Momentum().e()  / GeV; 
	double p_targ = 0.; 
	double e_targ = G4NucleiProperties::GetNuclearMass(targetNucleus.GetA_asInt(),targetNucleus.GetZ_asInt()) / GeV;
	double e_initial = e_proj + e_targ;
	// ... bug fix (March 2, 2020 - momentum per nucleon!)
	double a_proj = (double)(aTrack.GetDefinition()->GetAtomicMass()); // GetAtomicNumber());
	if(a_proj<1.0) a_proj = 1.0; // explanation: if particle is not an ion/proton, the GetAtomicMass returns 0
	double a_targ = (double)(targetNucleus.GetA_asInt());
	


	//* DEBUG messages
	if(fPrintDebug){
		std::cout<<"\n\n\n\n\n\n\n=============================================="<<std::endl;
		std::cout<<"Start interaction"<<std::endl;
		std::cout<<"id_proj="<<id_proj<<std::endl;
		std::cout<<"id_targ="<<id_targ<<std::endl;
		std::cout<<"p_proj="<<p_proj<<std::endl;
		std::cout<<"p_targ="<<p_targ<<std::endl;
	}
	
	
	// set up input particle type and energy
	fInterface->crmc_set(
		1,                  //fNCollision,
		p_proj / a_proj,    //fCfg.fProjectileMomentum (per nucleon!!!),
		p_targ / a_targ,    //fCfg.fTargetMomentum (per nucleon!!!),
		id_proj,            //fCfg.fProjectileId,
		id_targ);           //fCfg.fTargetId);
	
	//=================================================
	// sample 1 interaction until the energy 
	// conservation is fulfilled
	int resample_attampts = 1;
	double max_energy_diff = ENERGY_NON_CONSERVATION_EMAX_GEV;
	double energy_diff_coef =1.;
	double forbid_energy_corr = false;
	while(true){
		// run one interaction
		fInterface->crmc_generate(
			fTypeOutput, // fCfg.fTypoaut,
			1,           // iColl+1,
			gCRMC_data.fNParticles,
			gCRMC_data.fImpactParameter,
			gCRMC_data.fPartId[0],
			gCRMC_data.fPartPx[0],
			gCRMC_data.fPartPy[0],
			gCRMC_data.fPartPz[0],
			gCRMC_data.fPartEnergy[0],
			gCRMC_data.fPartMass[0],
			gCRMC_data.fPartStatus[0]);

		// split Z=0 A>1 "particles" into multiple neutrons
		SplitMultiNeutrons(gCRMC_data);

		// energy check
		double e_final =0;
		for(int i=0; i<gCRMC_data.fNParticles;i++){
			if (gCRMC_data.fPartStatus[i]!=1) continue; // only final state particles
			G4ParticleDefinition* pdef;
			int Z_test = (gCRMC_data.fPartId[i] - CRMC_ION_COEF_0)/CRMC_ION_COEF_Z;
			int A_test = (gCRMC_data.fPartId[i] - CRMC_ION_COEF_0 - CRMC_ION_COEF_Z*Z_test)/CRMC_ION_COEF_A;
			if(fPrintDebug){
				std::cout<<std::endl;
				std::cout<<"**********************************************************************"<<std::endl;
				std::cout<<"PDG test: " << gCRMC_data.fPartId[i] << std::endl;
				std::cout<<"fIonTable->GetIon(Z_test, A_test)                  = " << fIonTable->GetIon(Z_test, A_test)  << std::endl;
				std::cout<<"ParticleTable->FindParticle(gCRMC_data.fPartId[i]) = " << fParticleTable->FindParticle(gCRMC_data.fPartId[i])  << std::endl;
				std::cout<<"**********************************************************************"<<std::endl;
			}

			//pdef = fParticleTable->FindParticle(gCRMC_data.fPartId[i]);
			int pdef_errorcode;
			pdef = GetParticleDefinition(gCRMC_data.fPartId[i],pdef_errorcode);
			if(!pdef && IGNORE_PARTICLE_UNKNOWN_PDGID){
				continue;
			}

			double p2 = std::pow(gCRMC_data.fPartPx[i],2) + std::pow(gCRMC_data.fPartPy[i],2) +  std::pow(gCRMC_data.fPartPz[i],2);
			double mass  = pdef->GetPDGMass()/GeV;
			e_final  += std::sqrt(mass*mass + p2);
		}

		// Check if we need to resample again...
		double diff = fabs(e_final - e_initial);
		if(e_final!=0. && e_initial!=0. && USE_ENERGY_CORR) energy_diff_coef = e_final / e_initial;
		if(fPrintDebug){
			std::cout<< "# e_initial = " << e_initial << " GeV" << std::endl;
			std::cout<< "# e_final   = " << e_final   << " GeV" << std::endl;
			std::cout<< "# energy_diff_coef = " << energy_diff_coef << std::endl;
		}

		// energy conservation check, if yes
		if(!ENERGY_NON_CONSERVATION_RESAMPLE){
			// ===== NOCHECK ========== NOCHECK ============== NOCHECK ========
			break;
			// ===== NOCHECK ========== NOCHECK ============== NOCHECK ========
		}
		else if(diff<max_energy_diff || diff/e_initial < ENERGY_NON_CONSERVATION_FRACTION_MAX){
			// ===== OK ========== OK ============== OK ========
			forbid_energy_corr = true;
			break; // everything is fine, no need to resample, break the re-sampling loop
			// ===== OK ========== OK ============== OK ========
		}
		else if (resample_attampts<ENERGY_NON_CONSERVATION_FRACTION_MAX_ATTEMPT){
			resample_attampts++;
			std::cout<< std::endl;
			std::cout<< "#==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ====#" << std::endl;
			std::cout<< "#                                                                                  #" << std::endl;
			std::cout<< "# [HadronicInelasticModelCRMC::ApplyYourself]: Energy non conservation detected: #" << std::endl;
			std::cout<< "# e_initial = " << e_initial << " GeV"                                                << std::endl;
			std::cout<< "# e_final   = " << e_final   << " GeV"                                                << std::endl;
			std::cout<< "# diff      = " << diff      << " GeV"                                                << std::endl;
			std::cout<< "# Running attempt #" << resample_attampts                                             << std::endl;
			std::cout<< "#                                                                                  #" << std::endl;
			std::cout<< "#==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ====#" << std::endl;
			std::cout<< std::endl;
		}
		else if (max_energy_diff<ENERGY_NON_CONSERVATION_FRACTION_MAX_ENERGYTRY){
			std::cout<< std::endl;
			std::cout<< "#==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ====#" << std::endl;
			std::cout<< "# reached maximum number of attempts = " << ENERGY_NON_CONSERVATION_FRACTION_MAX_ATTEMPT << " ==> increasing twice the energy threshold!" << std::endl;
			max_energy_diff *= 2.;
			std::cout<< "# max_energy_diff = " << max_energy_diff << std::endl;
			std::cout<< "#==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ====#" << std::endl;
			std::cout<< std::endl;
			resample_attampts = 1;
		}
		else{
			std::cout<< std::endl;
			std::cout<< "#==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ====#" << std::endl;
			std::cout<< "# reached maximum number of attempts = " << ENERGY_NON_CONSERVATION_FRACTION_MAX_ATTEMPT << "not resampling any more!"  << std::endl;
			std::cout<< "#==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ==== WARNING ====#" << std::endl;
			std::cout<< std::endl;
			// ===== FAIL ========== FAIL ============== FAIl ========
			break;
			// ===== FAIL ========== FAIL ============== FAIl ========
		}
	}
	// ... finished sampling one interaction
	//=================================================
	
	// ... for DEBUG messages
	double totalenergy = 0;
	double totalz      = 0;
	double eaftertest0 = 0.;
	double eaftertest1 = 0.;
	double eaftertest2 = 0.;

	// save secondary particles for outputa
	for(int i=0; i<gCRMC_data.fNParticles;i++){ 
		//* Keep only final state particles
		// .. (-9 is the beam, 2 is a particle which decayed and 1 is final)
		if (gCRMC_data.fPartStatus[i]!=1) continue;
		
		if(fPrintDebug){
			std::cout<<"\n\nSecondary:"<< std::endl <<
				gCRMC_data.fPartId[i] << std::endl <<
				gCRMC_data.fPartPx[i] << std::endl <<
				gCRMC_data.fPartPy[i] << std::endl <<
				gCRMC_data.fPartPz[i] << std::endl <<
				gCRMC_data.fPartEnergy[i] << "  ENERGY  " << std::endl;
		}

		//G4ParticleDefinition* pdef = fParticleTable->FindParticle(gCRMC_data.fPartId[i]);
		int pdef_errorcode;
		G4ParticleDefinition* pdef = GetParticleDefinition(gCRMC_data.fPartId[i],pdef_errorcode);
		if(!pdef){
		  if(IGNORE_PARTICLE_UNKNOWN_PDGID){
			std::cout<<std::endl;
			std::cout<<"********************************************************************************************************"<<std::endl;
			std::cout<<" -- WARNING ----------------------------------------------------------------------------------- WARNING --  "<<std::endl;
			std::cout<<" [HadronicInelasticModelCRMC] Geant4 could not find particle definition for PDG ID = " << gCRMC_data.fPartId[i] << std::endl; 
			std::cout<<" [HadronicInelasticModelCRMC] Ignoring this particle. This might cause energy non-conservation!"      << std::endl; 
			std::cout<<" -- WARNING ----------------------------------------------------------------------------------- WARNING --  "<<std::endl;
			std::cout<<"********************************************************************************************************"<<std::endl;
			continue;
		  }else{
			std::cout<<std::endl;
			std::cout<<"********************************************************************************************************"<<std::endl;
			std::cout<<" -- ERROR ----------------------------------------------------------------------------------- ERROR --  "<<std::endl;
			std::cout<<" [HadronicInelasticModelCRMC] Geant4 could not find particle definition for PDG ID = " << gCRMC_data.fPartId[i] << std::endl; 
			std::cout<<" [HadronicInelasticModelCRMC] Throwing exception! Please report to: " << ERROR_REPORT_EMAIL   << std::endl; 
			std::cout<<" -- ERROR ----------------------------------------------------------------------------------- ERROR --  "<<std::endl;
			std::cout<<"********************************************************************************************************"<<std::endl;
			throw;
		  }
		}

		double part_e_corr = 1.;
		double part_p_corr = 1.;
		if(USE_ENERGY_CORR && !forbid_energy_corr && energy_diff_coef!=0){
			part_e_corr = 1./energy_diff_coef;
			double pbefore2 = std::pow(gCRMC_data.fPartPx[i],2) + std::pow(gCRMC_data.fPartPy[i],2) +  std::pow(gCRMC_data.fPartPz[i],2);
			double mass2    = std::pow(pdef->GetPDGMass()/GeV,2);  //std::pow(gCRMC_data.fPartEnergy[i],2) - pbefore2;
			double ebefore2 = pbefore2 + mass2;
			double pafter2  = ebefore2 * part_e_corr * part_e_corr - mass2;
			if(pbefore2) part_p_corr = std::sqrt(std::fabs(pafter2/pbefore2));
			if(fPrintDebug) std::cout<< "part_p_corr="<< part_p_corr << std::endl;
			eaftertest0 += std::sqrt(mass2 + pbefore2);
			eaftertest1 += std::sqrt(mass2 + pafter2);
		}
		
		G4DynamicParticle* part    = new G4DynamicParticle(pdef,G4ThreeVector(
			gCRMC_data.fPartPx[i]*GeV * part_p_corr,
			gCRMC_data.fPartPy[i]*GeV * part_p_corr,
			gCRMC_data.fPartPz[i]*GeV * part_p_corr
		));
		eaftertest2 += part->GetTotalEnergy ();
		finalState->AddSecondary(part);
		totalenergy += gCRMC_data.fPartEnergy[i];
		totalz += gCRMC_data.fPartPz[i];

	}

	if(fPrintDebug){
		std::cout<<"totalenergy (GeV) = " << totalenergy<<std::endl;
		std::cout<<"totalz (GeV)      = " << totalz<<std::endl;
		std::cout<<"initialz (GeV)    = " << p_proj + p_targ <<std::endl;
		std::cout<<"eaftertest0       = " << eaftertest0 <<std::endl;
		std::cout<<"eaftertest1       = " << eaftertest1 <<std::endl;
		std::cout<<"eaftertest2       = " << eaftertest2 <<std::endl;
		std::cout<<"Finishing interaction: "<<std::endl;
		const G4LorentzVector & p1 = aTrack.Get4Momentum ();
		std::cout<< "e=" << p1.e()<< " px=" <<p1.px() << " py=" << p1.py() << " pz="<<p1.pz() << std::endl;
		std::cout<< aTrack.GetDefinition()->GetAtomicNumber() << std::endl;
		std::cout<< aTrack.GetDefinition()->GetPDGCharge() << std::endl;
		std::cout<< targetNucleus.GetA_asInt() << std::endl;
		std::cout<< targetNucleus.GetZ_asInt() << std::endl;
		std::cout<<"Stop interaction"<<std::endl;
		std::cout<<"==============================================\n\n\n\n\n\n"<<std::endl;
	}
	//std::cout<<"finalState->GetNumberOfSecondaries()="<<finalState->GetNumberOfSecondaries()<< std::endl; // Debugging info
	return finalState;
}

G4bool HadronicInelasticModelCRMC::IsApplicable (const G4HadProjectile &, G4Nucleus &){
	return true;
}

G4ParticleDefinition* HadronicInelasticModelCRMC::GetParticleDefinition(long particle_id,int& error_code){ 
	G4ParticleDefinition* pdef = fParticleTable->FindParticle(particle_id);
	if(!pdef && particle_id > CRMC_ION_COEF_0){
		int Z = (particle_id - CRMC_ION_COEF_0)/CRMC_ION_COEF_Z;
		int A = (particle_id - CRMC_ION_COEF_0 - CRMC_ION_COEF_Z*Z)/CRMC_ION_COEF_A;
		if(IsMultiNeutron(Z,A)){
			error_code = PARTICLE_MULTI_NEUTRONS_ERRORCODE;
			pdef = NULL;
		}
		else{
			pdef = fIonTable->GetIon(Z, A);
		}
	}
	return pdef;
}

bool HadronicInelasticModelCRMC::IsMultiNeutron(int Z, int A){
	bool result = false;
	if (!Z && A>1){ 
		if(A<= SPLIT_MULTI_NEUTRONS_MAXN){
			result = true;
		}
		else{
			std::cout<<" [HadronicInelasticModelCRMC::IsMultiNeutron] ERROR A="<<
				A<<" is higher than "<< SPLIT_MULTI_NEUTRONS_MAXN << 
				" throwing exception!" << std::endl;
			throw;
		}
	}
	return result;
}

void HadronicInelasticModelCRMC::SplitMultiNeutrons(CRMCdata& CRMC_data){
	for(int i=0; i<CRMC_data.fNParticles;i++){
		// check if it is a final-state secondary particle
		if (CRMC_data.fPartStatus[i]!=1) continue; 

		int pdef_errorcode;
		GetParticleDefinition(CRMC_data.fPartId[i],pdef_errorcode);
		if(pdef_errorcode!=PARTICLE_MULTI_NEUTRONS_ERRORCODE) continue;

		//
		int particle_id = gCRMC_data.fPartId[i];
		int Z = (particle_id - CRMC_ION_COEF_0)/CRMC_ION_COEF_Z;
		int A = (particle_id - CRMC_ION_COEF_0 - CRMC_ION_COEF_Z*Z)/CRMC_ION_COEF_A;
		if(Z!=0 || A<2){
			std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] ERROR consistency check failed! Throwing exception! " << std::endl;
			throw;
		}

		//
		std::cout<<std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO splitting the floowing particle into neutrons: " << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    Z = " << Z << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    A = " << A << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    fPartId = " << CRMC_data.fPartId[i] << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    fPartPx = " << CRMC_data.fPartPx[i] << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    fPartPy = " << CRMC_data.fPartPy[i] << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    fPartPz = " << CRMC_data.fPartPz[i]  << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    fPartEnergy = " << CRMC_data.fPartEnergy[i]  << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    fPartMass   = " << CRMC_data.fPartMass[i]  << std::endl;
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] INFO    fPartStatus = " << CRMC_data.fPartStatus[i]  << std::endl;
		
		//
		int NEUTRON_PDG_ID = 2112;
		G4ParticleDefinition* p_n_def = fParticleTable->FindParticle(NEUTRON_PDG_ID);
		double m_n   = p_n_def->GetPDGMass()/GeV;
		double e_n   = CRMC_data.fPartEnergy[i]/A;
		int status_n = CRMC_data.fPartStatus[i];
		if(e_n<m_n){
			std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] WARNING neutron energy " << 
				e_n << " lower than neutron mass " <<
				m_n << " assigning e_n = m_n     " << std::endl;
			e_n = m_n;
		}
		double p_tot_before = std::sqrt(
			CRMC_data.fPartPx[i]*CRMC_data.fPartPx[i] + 
			CRMC_data.fPartPy[i]*CRMC_data.fPartPy[i] + 
			CRMC_data.fPartPz[i]*CRMC_data.fPartPz[i] 
		);
		double p_tot_after = std::sqrt(e_n*e_n - m_n*m_n);
		double px_n = 0;
		double py_n = 0;
		double pz_n = 0;
		if (p_tot_before>0. && p_tot_after>0.){
			px_n = CRMC_data.fPartPx[i] * p_tot_after / p_tot_before;
			py_n = CRMC_data.fPartPy[i] * p_tot_after / p_tot_before;
			pz_n = CRMC_data.fPartPz[i] * p_tot_after / p_tot_before;
		}
		for(int j=0;j<A;j++){
			int i_neutron = j ? CRMC_data.fNParticles+j : i;
			CRMC_data.fPartId[i_neutron] = NEUTRON_PDG_ID;
			CRMC_data.fPartPx[i_neutron] = px_n;
			CRMC_data.fPartPy[i_neutron] = py_n;
			CRMC_data.fPartPz[i_neutron] = pz_n;
			CRMC_data.fPartEnergy[i_neutron] =  e_n;
			CRMC_data.fPartMass[i_neutron] = m_n;
			CRMC_data.fPartStatus[i_neutron] = status_n;
		}
		CRMC_data.fNParticles+=A-1;

		//
		std::cout<<" [HadronicInelasticModelCRMC::SplitMultiNeutrons] done for a particle. " << std::endl;
		std::cout<<std::endl;
	}
}

#endif //G4_USE_CRMC
