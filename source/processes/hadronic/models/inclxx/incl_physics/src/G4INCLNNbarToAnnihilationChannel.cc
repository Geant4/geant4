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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"
#include "G4EnvironmentUtils.hh"
#include "G4INCLNNbarToAnnihilationChannel.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"

namespace G4INCL {
	
	NNbarToAnnihilationChannel::NNbarToAnnihilationChannel(Nucleus *n, Particle *p1, Particle *p2)
		:theNucleus(n), particle1(p1), particle2(p2)
		{}
	
	NNbarToAnnihilationChannel::~NNbarToAnnihilationChannel(){}
	
//fill probabilities and particle types from datafile and return probability sum for normalization
G4double NNbarToAnnihilationChannel::read_file(std::string filename, std::vector<G4double>& probabilities, std::vector<std::vector<std::string>>& particle_types) {
  std::ifstream file(filename);
  G4double sum_probs = 0.0;
  if (file.is_open()) {
      std::string line;
      while (getline(file, line)) {
          std::istringstream iss(line);
          G4double prob;
          iss >> prob;
          sum_probs += prob;
          probabilities.push_back(prob);
          std::vector<std::string> types;
          std::string type;
          while (iss >> type) {
              types.push_back(type);
          }
          particle_types.push_back(types);
      }
  }
  else std::cout << "ERROR no fread_file " << filename << std::endl;
  return sum_probs;
}

//this function will tell you the FS line number from the datafile based on your random value
G4int NNbarToAnnihilationChannel::findStringNumber(G4double rdm, std::vector<G4double> yields) {
		G4int stringNumber = -1;
		G4double smallestsum = 0.0;
		G4double biggestsum = yields[0];
		//std::cout << "initial input " << rdm << std::endl;
		for (G4int i = 0; i < static_cast<G4int>(yields.size()-1); i++) {
		    if (rdm >= smallestsum && rdm <= biggestsum) {
		        //std::cout << smallestsum << " and " << biggestsum << std::endl;
		        stringNumber = i+1;
		    }
		    smallestsum += yields[i];
		    biggestsum += yields[i+1];
		}
    if(stringNumber==-1) stringNumber = static_cast<G4int>(yields.size());
		if(stringNumber==-1){
		  INCL_ERROR("ERROR in findStringNumber (stringNumber=-1)");
		  std::cout << "ERROR in findStringNumber" << std::endl;
		}
		return stringNumber;
}


void NNbarToAnnihilationChannel::fillFinalState(FinalState *fs) {

		Particle *nucleon;
		Particle *antinucleon;
		
		if(particle1->isNucleon()){
			nucleon = particle1;
			antinucleon = particle2;
		}
		else{
			nucleon = particle2;
			antinucleon = particle1;
		}

		const G4double plab = 0.001*KinematicsUtils::momentumInLab(particle1, particle2); //GeV
		const G4double sqrtS = KinematicsUtils::totalEnergyInCM(nucleon, antinucleon);
		G4double rdm = Random::shoot();

    const std::vector<G4double> BFMM6 = {66.098, 0.153, -4.576, -38.319, 6.625}; //ppbar annihilation xs
    const std::vector<G4double> BFMM1 = {119.066, 6.251, -0.006, -60.046, 11.958}; //ppbar total xs
    const std::vector<G4double> BFMM471 = {108.104, 15.708, 0.832, -54.632, -6.958}; //npbar total xs

  //PPbar annihilation xs
    const std::vector<G4double> PPbar_pip_pim = {0.637, -0.340, -0.003, -0.439, 0.144};
    const std::vector<G4double> PPbar_pip_pim_pi0 = {-2.065, 4.893, -1.130, 1.231, -0.212};
		const std::vector<G4double> PPbar_pip_pim_omega = {3.020, 0.425, -0.029, -3.420, 0.867};
		const std::vector<G4double> PPbar_pip_pim_Kp_Km = {-1.295, 1.897, -0.001, -0.365, 0.044};
		const std::vector<G4double> PPbar_pip_pim_pi0_Kp_Km = {-12.220, 12.509, -0.351, 4.682, -0.777};
		const std::vector<G4double> PPbar_2pip_2pim = {3.547, 0.095, 0.957, -3.444, 0.685};
		const std::vector<G4double> PPbar_2pip_2pim_pi0 = {13.044, 1.449, 0.695, -12.313, 1.627};
		const std::vector<G4double> PPbar_2pip_2pim_3pi0 = {6.398, 0.199, -1.103, -1.271, -0.380};
		const std::vector<G4double> PPbar_3pip_3pim = {1.490, 0.240, 0.002, -1.012, 0.134};
		const std::vector<G4double> PPbar_3pip_3pim_pi0 = {0.286, 1.634, -1.369, 3.099, -1.294};
		const std::vector<G4double> PPbar_3pip_3pim_2pi0 = {-11.370, 12.503, -0.680, 10.059, -2.501};
		const std::vector<G4double> PPbar_3pip_3pim_3pi0 = {-14.732, 12.338, -0.724, 11.342, -2.224};
		const std::vector<G4double> PPbar_4pip_4pim = {-1.574, 1.607, -0.864, 1.253, -0.276};
		const std::vector<G4double> PPbar_4pip_4pim_pi0 = {-1.096, 0.977, -0.995, 1.007, -0.171};

	//NPbar annihilation xs
    const std::vector<G4double> NPbar_pip_2pim = {-12.116, 14.485, -0.094, -1.632, 0.882, 5.000};
		const std::vector<G4double> NPbar_pip_2pim_2pi0 = {8.276, 5.057, 0.483, -15.864, 2.552, 7.000};
		const std::vector<G4double> NPbar_pip_2pim_3pi0 = {-1.500, 9.574, 0.528, -11.633, -0.615, 7.000};
		const std::vector<G4double> NPbar_pip_2pim_pi0 = {7.999, 4.135, 0.608, -14.136, 1.590, 7.000};
		const std::vector<G4double> NPbar_pip_pim_pi0_Km_K0 = {0.083, 0.091, -1.709, 0.284, -0.107};
		const std::vector<G4double> NPbar_pip_pim_Km_K0 = {0.003, 0.297, -0.001, -0.143, 0.052};
		const std::vector<G4double> NPbar_2pip_3pim_pi0 = {-14.701, 22.258, -0.001, -3.094, -0.190};
		const std::vector<G4double> NPbar_2pip_3pim = {-0.616, 4.575, -0.002, -1.921, -0.153};


		// File names
		    #ifdef INCLXX_IN_GEANT4_MODE
		       if(!G4FindDataDir("G4INCLDATA")) {
		        G4ExceptionDescription ed;
		        ed << " Data missing: set environment variable G4INCLDATA\n"
		           << " to point to the directory containing data files needed\n"
		           << " by the INCL++ model" << G4endl;
		           G4Exception("G4INCLDataFile::readData()","inflightppbarFS.dat, ...",
		                FatalException, ed);
		      }
		      G4String dataPath0{G4FindDataDir("G4INCLDATA")};
		      G4String dataPathppbar(dataPath0 + "/inflightppbarFS.dat");
		      G4String dataPathnpbar(dataPath0 + "/inflightnpbarFS.dat");
		      G4String dataPathppbark(dataPath0 + "/inflightppbarFSkaonic.dat");
		      G4String dataPathnpbark(dataPath0 + "/inflightnpbarFSkaonic.dat");
		      G4String dataPathpnbar(dataPath0 + "/inflightpnbarFS.dat"); //nbar case
		      G4String dataPathpnbark(dataPath0 + "/inflightpnbarFSkaonic.dat"); // nbar case
		    #else
		      //Config *theConfig = new G4INCL::Config;
			  //theConfig->setINCLXXDataFilePath(G4INCL::theINCLXXDataFilePath);
		      Config const *theConfig=theNucleus->getStore()->getConfig();
		      std::string path;
		      if(theConfig)
		        path = theConfig->getINCLXXDataFilePath();
		      std::string dataPathppbar(path + "/inflightppbarFS.dat");
		      INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N ppbar final states" << dataPathppbar << '\n');
		      std::string dataPathnpbar(path + "/inflightnpbarFS.dat");
		      INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N npbar final states" << dataPathnpbar << '\n');
		      std::string dataPathppbark(path + "/inflightppbarFSkaonic.dat");
		      INCL_DEBUG("Reading https://doi.org/10.1016/j.physrep.2005.03.002 ppbar kaonic final states" << dataPathppbark << '\n');
		      std::string dataPathnpbark(path + "/inflightnpbarFSkaonic.dat");
		      INCL_DEBUG("Reading https://doi.org/10.1007/BF02818764 and https://link.springer.com/article/10.1007/BF02754930 npbar kaonic final states" << dataPathnpbark << '\n');
		      std::string dataPathpnbar(path + "/inflightpnbarFS.dat"); //  nbar case
		      std::string dataPathpnbark(path + "/inflightpnbarFSkaonic.dat"); //  nbar case
		   #endif
		    /*std::string path = {"/home/zdemid/INCL/inclcode/data"};
		    std::string dataPathppbar(path + "/inflightppbarFS.dat");
		    INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N ppbar final states" << dataPathppbar << '\n');
		    std::string dataPathnpbar(path + "/inflightnpbarFS.dat");
		    INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N npbar final states" << dataPathnpbar << '\n');
		    std::string dataPathppbark(path + "/inflightppbarFSkaonic.dat");
		    INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N ppbar kaonic final states" << dataPathppbark << '\n');
		    std::string dataPathnpbark(path + "/inflightnpbarFSkaonic.dat");
		    INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N npbar kaonic final states" << dataPathnpbark << '\n');
		    every time we remove lines for which we have data from BFMM
		    */

		std::vector<G4double> probabilities; //will store each FS yield
    std::vector<std::vector<std::string>> particle_types; //will store particle names
    G4double sum; //will contain a sum of probabilities of all FS in the file
		const G4double kaonicFSprob=0.05; //probability to kave kaonic FS


    ParticleList list;
    //list.push_back(nucleon);
		//list.push_back(antinucleon);
		// NNbar will not be in the list because they annihilate
    const ThreeVector &rcol = nucleon->getPosition();
		const ThreeVector zero;
		
		//setting types of new particles and pushing them back to the list
		if(nucleon->getType()==Neutron && antinucleon->getType()==antiProton){
			//std::cout << "npbar"<< std::endl;
			const G4double totalpnbar = KinematicsUtils::compute_xs(BFMM6, plab)*KinematicsUtils::compute_xs(BFMM471, plab)/KinematicsUtils::compute_xs(BFMM1, plab);
			// xs is same for npbar, but the fs has different charge

			if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab)) {
			    // First condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiMinus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab)) {
			    // Second condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(PiZero, zero, rcol);
			    Particle* p5 = new Particle(PiMinus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab)) {
			    // Third condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(PiZero, zero, rcol);
			    Particle* p5 = new Particle(PiZero, zero, rcol);
			    Particle* p6 = new Particle(PiMinus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab)) {
			    // Fourth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_pi0_Km_K0, plab)) {
			    // Fifth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(KMinus, zero, rcol);
			    Particle* p5 = new Particle(KZero, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_pi0_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_Km_K0, plab)) {
			    // Sixth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(KMinus, zero, rcol);
			    Particle* p4 = new Particle(KZero, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_pi0_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_2pip_3pim_pi0, plab)) {
			    // Seventh condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiPlus, zero, rcol);
			    Particle* p3 = new Particle(PiMinus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiMinus, zero, rcol);
			    Particle* p6 = new Particle(PiZero, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_pi0_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_2pip_3pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_2pip_3pim, plab)) {
			    // Eighth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiPlus, zero, rcol);
			    Particle* p3 = new Particle(PiMinus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiMinus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			} else {
			    // Default condition    
				//std::cout << "default condition pnbar"<< std::endl;
				if(rdm < (1.-kaonicFSprob)){ // pionic/kaonic choice
			        INCL_DEBUG("pionic npbar final state chosen" << '\n');
			        sum = read_file(dataPathnpbar, probabilities, particle_types);
			        rdm = (rdm/(1.-kaonicFSprob))*sum; //normalize by the sum of probabilities in the file
			        //now get the line number in the file where the FS particles are stored:
			        G4int n = findStringNumber(rdm, probabilities)-1;
			        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
			          if(particle_types[n][j] == "pi0"){
			            Particle *p = new Particle(PiZero, zero, rcol);
			            list.push_back(p);
			          }
			          else if(particle_types[n][j] == "pi-"){
			            Particle *p = new Particle(PiMinus, zero, rcol);
			            list.push_back(p);
			          }
			          else if(particle_types[n][j] == "pi+"){
			            Particle *p = new Particle(PiPlus, zero, rcol);
			            list.push_back(p);
			          }
			          else if(particle_types[n][j] == "omega"){
			            Particle *p = new Particle(Omega, zero, rcol);
			            list.push_back(p);
			          }
			          else if(particle_types[n][j] == "eta"){
			            Particle *p = new Particle(Eta, zero, rcol);
			            list.push_back(p);
			          }
			          else{
			            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
			            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
			              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
			            }
			          }
        			}
      			} // end of pionic option
		        else{
			        INCL_DEBUG("kaonic npbar final state chosen" << '\n');
			        sum = read_file(dataPathnpbark, probabilities, particle_types);
			        rdm = ((1-rdm)/kaonicFSprob)*sum;//3837 normalize by the sum of probabilities in the file
			        //now get the line number in the file where the FS particles are stored:
			        G4int n = findStringNumber(rdm, probabilities)-1;
			        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
				        if(particle_types[n][j] == "pi0"){
				            Particle *p = new Particle(PiZero, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "pi-"){
				            Particle *p = new Particle(PiMinus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "pi+"){
				            Particle *p = new Particle(PiPlus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "omega"){
				            Particle *p = new Particle(Omega, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "eta"){
				            Particle *p = new Particle(Eta, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K-"){
				            Particle *p = new Particle(KMinus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K+"){
				            Particle *p = new Particle(KPlus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K0"){
				            Particle *p = new Particle(KZero, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K0b"){
				            Particle *p = new Particle(KZeroBar, zero, rcol);
				            list.push_back(p);
				          }
				        else{
				            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
				            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
				              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
				            }
				        }
		    		}
    			} // end of kaonic option
			} // end of default annihilation

		}
		else if(nucleon->getType()==Proton && antinucleon->getType()==antiNeutron){
			const G4double totalpnbar = KinematicsUtils::compute_xs(BFMM6, plab)*KinematicsUtils::compute_xs(BFMM471, plab)/KinematicsUtils::compute_xs(BFMM1, plab);
			// xs is same for npbar, but the fs has different charge

			if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab)) {
			    // First condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab)) {
			    // Second condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(PiZero, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab)) {
			    // Third condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(PiZero, zero, rcol);
			    Particle* p5 = new Particle(PiZero, zero, rcol);
			    Particle* p6 = new Particle(PiPlus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab)) {
			    // Fourth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(PiPlus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_pi0_Km_K0, plab)) {
			    // Fifth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(KPlus, zero, rcol);
			    Particle* p5 = new Particle(KZeroBar, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_pi0_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_Km_K0, plab)) {
			    // Sixth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(KPlus, zero, rcol);
			    Particle* p4 = new Particle(KZeroBar, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_pi0_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_2pip_3pim_pi0, plab)) {
			    // Seventh condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiPlus, zero, rcol);
			    Particle* p3 = new Particle(PiMinus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    Particle* p6 = new Particle(PiZero, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			} else if (rdm * totalpnbar < KinematicsUtils::compute_xs(NPbar_pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_pi0_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_pip_pim_Km_K0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_2pip_3pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(NPbar_2pip_3pim, plab)) {
			    // Eighth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiPlus, zero, rcol);
			    Particle* p3 = new Particle(PiMinus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    
			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			} else {
			    // Default condition    
				if(rdm < (1.-kaonicFSprob)){ // pionic/kaonic choice
			        INCL_DEBUG("pionic pnbar final state chosen" << '\n');
			        sum = read_file(dataPathpnbar, probabilities, particle_types);
			        rdm = (rdm/(1.-kaonicFSprob))*sum; //99.95 normalize by the sum of probabilities in the file
			        //now get the line number in the file where the FS particles are stored:
			        G4int n = findStringNumber(rdm, probabilities)-1;
			        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
			          if(particle_types[n][j] == "pi0"){
			            Particle *p = new Particle(PiZero, zero, rcol);
			            list.push_back(p);
			          }
			          else if(particle_types[n][j] == "pi-"){
			            Particle *p = new Particle(PiMinus, zero, rcol);
			            list.push_back(p);
			          }
			          else if(particle_types[n][j] == "pi+"){
			            Particle *p = new Particle(PiPlus, zero, rcol);
			            list.push_back(p);
			          }
			          else if(particle_types[n][j] == "omega"){
			            Particle *p = new Particle(Omega, zero, rcol);
			            list.push_back(p);
			          }
			          else if(particle_types[n][j] == "eta"){
			            Particle *p = new Particle(Eta, zero, rcol);
			            list.push_back(p);
			          }
			          else{
			            INCL_ERROR("Some non-existing FS particle detected when reading nbar FS files");
			            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
			              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
			            }
			          }
        			}
      			} // end of pionic option
		        else{
			        INCL_DEBUG("kaonic pnbar final state chosen" << '\n');
			        sum = read_file(dataPathnpbark, probabilities, particle_types);
			        rdm = ((1-rdm)/kaonicFSprob)*sum;//3837 normalize by the sum of probabilities in the file
			        //now get the line number in the file where the FS particles are stored:
			        G4int n = findStringNumber(rdm, probabilities)-1;
			        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
				        if(particle_types[n][j] == "pi0"){
				            Particle *p = new Particle(PiZero, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "pi-"){
				            Particle *p = new Particle(PiMinus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "pi+"){
				            Particle *p = new Particle(PiPlus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "omega"){
				            Particle *p = new Particle(Omega, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "eta"){
				            Particle *p = new Particle(Eta, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K-"){
				            Particle *p = new Particle(KMinus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K+"){
				            Particle *p = new Particle(KPlus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K0"){
				            Particle *p = new Particle(KZero, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K0b"){
				            Particle *p = new Particle(KZeroBar, zero, rcol);
				            list.push_back(p);
				          }
				        else{
				            INCL_ERROR("Some non-existing FS particle detected when reading pnbar FS files");
				            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
				              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;      }
				        }
		    		}
    			} // end of kaonic option
			} // end of default annihilation

		}
		else{ //ppbar or nnbar
			//std::cout << "ppbar or nnbar"<< std::endl;
			const G4double totalppbar = KinematicsUtils::compute_xs(BFMM6, plab);
			// same for nnbar

			if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab)) {
			    // First condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			
			
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab)) {
			    // Second condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab)) {
			    // Third condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(Omega, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab)) {
			    // Fourth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(KPlus, zero, rcol);
			    Particle* p4 = new Particle(KMinus, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab)) {
			    // Fifth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiZero, zero, rcol);
			    Particle* p4 = new Particle(KPlus, zero, rcol);
			    Particle* p5 = new Particle(KMinus, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab)) {
			    // Sixth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_pi0, plab)) {
			    // Seventh condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiZero, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab) +
                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_pi0, plab) +
                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_3pi0, plab)) {
			    // Eighth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiZero, zero, rcol);
			    Particle* p6 = new Particle(PiZero, zero, rcol);
			    Particle* p7 = new Particle(PiZero, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			    list.push_back(p7);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim, plab)) {
			    // Ninth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    Particle* p6 = new Particle(PiMinus, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_pi0, plab)) {
			    // Tenth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    Particle* p6 = new Particle(PiMinus, zero, rcol);
			    Particle* p7 = new Particle(PiZero, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			    list.push_back(p7);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_2pi0, plab)) {
			    // Eleventh condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    Particle* p6 = new Particle(PiMinus, zero, rcol);
			    Particle* p7 = new Particle(PiZero, zero, rcol);
			    Particle* p8 = new Particle(PiZero, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			    list.push_back(p7);
			    list.push_back(p8);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_3pi0, plab)) {
			    // Twelfth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    Particle* p6 = new Particle(PiMinus, zero, rcol);
			    Particle* p7 = new Particle(PiZero, zero, rcol);
			    Particle* p8 = new Particle(PiZero, zero, rcol);
			    Particle* p9 = new Particle(PiZero, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			    list.push_back(p7);
			    list.push_back(p8);
			    list.push_back(p9);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_4pip_4pim, plab)) {
			    // Thirteenth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    Particle* p6 = new Particle(PiMinus, zero, rcol);
			    Particle* p7 = new Particle(PiPlus, zero, rcol);
			    Particle* p8 = new Particle(PiMinus, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			    list.push_back(p7);
			    list.push_back(p8);
			} else if (rdm * totalppbar < KinematicsUtils::compute_xs(PPbar_pip_pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_omega, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_pip_pim_pi0_Kp_Km, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_2pip_2pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_2pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_3pip_3pim_3pi0, plab) +
			                            KinematicsUtils::compute_xs(PPbar_4pip_4pim, plab) +
			                            KinematicsUtils::compute_xs(PPbar_4pip_4pim_pi0, plab)) {
			    // Fourteenth condition
			    Particle* p1 = new Particle(PiPlus, zero, rcol);
			    Particle* p2 = new Particle(PiMinus, zero, rcol);
			    Particle* p3 = new Particle(PiPlus, zero, rcol);
			    Particle* p4 = new Particle(PiMinus, zero, rcol);
			    Particle* p5 = new Particle(PiPlus, zero, rcol);
			    Particle* p6 = new Particle(PiMinus, zero, rcol);
			    Particle* p7 = new Particle(PiPlus, zero, rcol);
			    Particle* p8 = new Particle(PiMinus, zero, rcol);
			    Particle* p9 = new Particle(PiZero, zero, rcol);

			    list.push_back(p1);
			    list.push_back(p2);
			    list.push_back(p3);
			    list.push_back(p4);
			    list.push_back(p5);
			    list.push_back(p6);
			    list.push_back(p7);
			    list.push_back(p8);
			    list.push_back(p9);
			} else {
			    // Default condition
			    if(rdm < (1.-kaonicFSprob)){ // pionic FS was chosen
			        INCL_DEBUG("pionic pp final state chosen" << '\n');
			        sum = read_file(dataPathppbar, probabilities, particle_types);
			        rdm = (rdm/(1.-kaonicFSprob))*sum; //99.88 normalize by the sum of probabilities in the file
			        //now get the line number in the file where the FS particles are stored:
			        G4int n = findStringNumber(rdm, probabilities)-1;
			        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
				        if(particle_types[n][j] == "pi0"){
				            Particle *p = new Particle(PiZero, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "pi-"){
				            Particle *p = new Particle(PiMinus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "pi+"){
				            Particle *p = new Particle(PiPlus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "omega"){
				            Particle *p = new Particle(Omega, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "eta"){
				            Particle *p = new Particle(Eta, zero, rcol);
				            list.push_back(p);
				          }
				        else{
							INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
							for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
								std::cout << "gotcha! " << particle_types[n][jj] << std::endl;      }
		        		  }
		    			}
				} //end of pionic option
				else{
			        INCL_DEBUG("kaonic pp final state chosen" << '\n');
			        sum = read_file(dataPathppbark, probabilities, particle_types);
			        rdm = ((1-rdm)/kaonicFSprob)*sum;//2670 normalize by the sum of probabilities in the file
			        //now get the line number in the file where the FS particles are stored:
			        G4int n = findStringNumber(rdm, probabilities)-1;
				    for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
				        if(particle_types[n][j] == "pi0"){
				            Particle *p = new Particle(PiZero, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "pi-"){
				            Particle *p = new Particle(PiMinus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "pi+"){
				            Particle *p = new Particle(PiPlus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "omega"){
				            Particle *p = new Particle(Omega, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "eta"){
				            Particle *p = new Particle(Eta, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K-"){
				            Particle *p = new Particle(KMinus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K+"){
				            Particle *p = new Particle(KPlus, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K0"){
				            Particle *p = new Particle(KZero, zero, rcol);
				            list.push_back(p);
				          }
				        else if(particle_types[n][j] == "K0b"){
				            Particle *p = new Particle(KZeroBar, zero, rcol);
				            list.push_back(p);
				          }
				        else{
				            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
				            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
				            	std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
				            }
				          }
				        }
      			} // end of kaonic option
					} // end of default condition
				} // end of ppbar and nnbar case



		nucleon->setType(list[0]->getType());
		antinucleon->setType(list[1]->getType());

		ParticleList finallist;

		finallist.push_back(nucleon);
		finallist.push_back(antinucleon);

		if(list.size() > 2){
			for (G4int i = 2; i < (G4int)(list.size()); i++) {
    		finallist.push_back(list[i]);
			}
		}

		if(finallist.size()==2){
			G4double mn=nucleon->getMass();
			G4double my=antinucleon->getMass();
			
			G4double ey=(sqrtS*sqrtS+my*my-mn*mn)/(2*sqrtS);
			G4double en=std::sqrt(ey*ey-my*my+mn*mn);
			nucleon->setEnergy(en);
			antinucleon->setEnergy(ey);
			G4double py=std::sqrt(ey*ey-my*my);

			ThreeVector mom_antinucleon = Random::normVector(py);
			antinucleon->setMomentum(mom_antinucleon);
			nucleon->setMomentum(-mom_antinucleon);
		}
		else if(finallist.size() > 2){
			PhaseSpaceGenerator::generate(sqrtS, finallist);
		}
		else{
			INCL_ERROR("less than 2 mesons in NNbar annihilation!" << '\n');
		}


		for (G4int i = 0; i < 2; i++) {
    	fs->addModifiedParticle(finallist[i]);
		}
		if(finallist.size()>2){
			for (G4int i = 2; i < (G4int)(list.size()); i++) {
    		fs->addCreatedParticle(finallist[i]);
			}
		}


		//fs->addDestroyedParticle(nucleon);
		//fs->addDestroyedParticle(antinucleon);

	
	}
}

