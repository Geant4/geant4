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

#include "G4INCLPbarAtrestEntryChannel.hh"
#include "G4INCLRootFinder.hh"
#include "G4INCLIntersection.hh"
#include "G4INCLCascade.hh"
#include <algorithm>
#include "G4INCLParticle.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"
#include <algorithm>
#include "G4INCLPhaseSpaceGenerator.hh"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "G4INCLHFB.hh"
#include "G4INCLParticleEntryAvatar.hh"
#include "G4INCLNuclearDensityFactory.hh"
#include "G4INCLNDFWoodsSaxon.hh"
#include "G4INCLNDFModifiedHarmonicOscillator.hh"
#include "G4INCLNDFGaussian.hh"
#include "G4INCLNDFParis.hh"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

namespace G4INCL {

  PbarAtrestEntryChannel::PbarAtrestEntryChannel(Nucleus *n, Particle *p)
    :theNucleus(n), theParticle(p)
  {}

  PbarAtrestEntryChannel::~PbarAtrestEntryChannel()
  {}
 
  //fill probabilities and particle types from datafile and return probability sum for normalization
  G4double PbarAtrestEntryChannel::read_file(std::string filename, std::vector<G4double>& probabilities, std::vector<std::vector<std::string>>& particle_types) {
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
  G4int PbarAtrestEntryChannel::findStringNumber(G4double rdm, std::vector<G4double> yields) {
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

  G4double PbarAtrestEntryChannel::fctrl(const G4double arg1){  
    G4double factorial=1.0; 
    for(G4int k=1; k<=arg1; k++){    
      factorial *= k;
    }        
    return factorial;
  }
  G4double PbarAtrestEntryChannel::r1(const G4int n) {
    return std::pow(fctrl(2.0*n),-0.5); 
  }
  G4double PbarAtrestEntryChannel::r2(const G4int n) {
    G4int Z = theNucleus->getZ();
    return std::pow((Z/(14.4*n)), 1.5);
  }
  G4double PbarAtrestEntryChannel::r3(G4double x, const G4int n) {
    G4int Z = theNucleus->getZ();
    return std::pow((x)*Z/(n*14.4),(n-1)); //why x is vector here?
  }
  G4double PbarAtrestEntryChannel::r4(G4double x, const G4int n) {
    G4int Z = theNucleus->getZ();
    return std::exp((-x*Z)/(n*28.8));
  }

  
  G4double PbarAtrestEntryChannel::radial_wavefunction(G4double x, const G4int n){
    return  r1(n)*r2(n)*r3(x,n)*r4(x,n); //Radial wave function par=(n, Z) 
  }

/*
  G4double PbarAtrestEntryChannel::densityP(G4double *x, G4double *par){
    return  0.16800136/(1.0+std::exp((x[0]-par[2])/par[3])); //P nuclear density
*/
  G4double PbarAtrestEntryChannel::densityP(G4double x) {
        
        const G4bool isProton = ProtonIsTheVictim();
        G4int Z = theNucleus->getZ(); //was modified in Cascade.cc
        G4int A = theNucleus->getA(); //was modified in Cascade.cc
        A++; //restoration of original A value before annihilation
        if(isProton == true){Z++;} //restoration of original Z value before annihilation

        if(A > 19) {
          G4double radius = ParticleTable::getRadiusParameter(Proton, A, Z);
          G4double diffuseness = ParticleTable::getSurfaceDiffuseness(Proton, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Proton, A, Z);
          NuclearDensityFunctions::WoodsSaxon rDensityFunction(radius, maximumRadius, diffuseness);
          return ((x!=0.) ? rDensityFunction(x)/(x*x) : 1.);
        } else if(A <= 19 && A > 6) {
          G4double radius = ParticleTable::getRadiusParameter(Proton, A, Z);
          G4double diffuseness = ParticleTable::getSurfaceDiffuseness(Proton, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Proton, A, Z);
          NuclearDensityFunctions::ModifiedHarmonicOscillator rDensityFunction(radius, maximumRadius, diffuseness);
          return ((x!=0.) ? rDensityFunction(x)/(x*x) : 1.);
        } else if(A <= 6 && A > 2) { // Gaussian distribution for light nuclei
          G4double radius = ParticleTable::getRadiusParameter(Proton, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Proton, A, Z);
          NuclearDensityFunctions::Gaussian rDensityFunction(maximumRadius, Math::oneOverSqrtThree * radius);
          return ((x!=0.) ? rDensityFunction(x)/(x*x) : 1.);
        } else if(A == 2 && Z == 1) { // density from the Paris potential for deuterons
          NuclearDensityFunctions::ParisR rDensityFunction;
          return ((x!=0.) ? rDensityFunction(x)/(x*x) : 1.);
        } else {
          INCL_ERROR("No nuclear density function for target A = "
                << A << " Z = " << Z << '\n');
          return 0.0;
        }

}
/*  
  }
  G4double PbarAtrestEntryChannel::densityN(G4double *x, G4double *par){
    return  0.16800136/(1.0+std::exp((x[0]-par[2])/par[3])); //N nuclear density
  }
*/      
  G4double PbarAtrestEntryChannel::densityN(G4double x) {
        
        const G4bool isProton = ProtonIsTheVictim();
        G4int Z = theNucleus->getZ(); //was modified in Cascade.cc
        G4int A = theNucleus->getA(); //was modified in Cascade.cc
        A++; //restoration of original A value before annihilation
        if(isProton == true){Z++;} //restoration of original Z value before annihilation

        if(A > 19) {
          G4double radius = ParticleTable::getRadiusParameter(Neutron, A, Z);
          G4double diffuseness = ParticleTable::getSurfaceDiffuseness(Neutron, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Neutron, A, Z);
          NuclearDensityFunctions::WoodsSaxon rDensityFunction(radius, maximumRadius, diffuseness);
          return ((x!=0.) ? rDensityFunction(x)/(x*x) : 1.);
        } else if(A <= 19 && A > 6) {
          G4double radius = ParticleTable::getRadiusParameter(Neutron, A, Z);
          G4double diffuseness = ParticleTable::getSurfaceDiffuseness(Neutron, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Neutron, A, Z);
          NuclearDensityFunctions::ModifiedHarmonicOscillator rDensityFunction(radius, maximumRadius, diffuseness);
          return ((x!=0.) ? rDensityFunction(x)/(x*x) : 1.);
        } else if(A <= 6 && A > 2) { // Gaussian distribution for light nuclei
          G4double radius = ParticleTable::getRadiusParameter(Neutron, A, Z);
          G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Neutron, A, Z);
          NuclearDensityFunctions::Gaussian rDensityFunction(maximumRadius, Math::oneOverSqrtThree * radius);
          return ((x!=0.) ? rDensityFunction(x)/(x*x) : 1.);
        } else if(A == 2 && Z == 1) { // density from the Paris potential for deuterons
          NuclearDensityFunctions::ParisR rDensityFunction;
          return ((x!=0.) ? rDensityFunction(x)/(x*x) : 1.);
        } else {
          INCL_ERROR("No nuclear density function for target A = "
                << A << " Z = " << Z << '\n');
          return 0.0;
        }

}


  G4double PbarAtrestEntryChannel::overlapP(G4double &x, const G4int n){
    return  x*x*r1(n)*r2(n)*r3(x,n)*r4(x,n)*r1(n)*r2(n)*r3(x,n)*r4(x,n)*densityP(x);
  }   
  G4double PbarAtrestEntryChannel::overlapN(G4double &x, const G4int n){
    return  x*x*r1(n)*r2(n)*r3(x,n)*r4(x,n)*r1(n)*r2(n)*r3(x,n)*r4(x,n)*densityN(x);
  }
  

  ParticleList PbarAtrestEntryChannel::makeMesonStar() {//This function creates a set of mesons with momenta

    // File names
    #ifdef INCLXX_IN_GEANT4_MODE
       if(!G4FindDataDir("G4INCLDATA")) {
        G4ExceptionDescription ed;
        ed << " Data missing: set environment variable G4INCLDATA\n"
           << " to point to the directory containing data files needed\n"
           << " by the INCL++ model" << G4endl;
           G4Exception("G4INCLDataFile::readData()","rawppbarFS.dat, ...",
                FatalException, ed);
      }
      G4String dataPath0{G4FindDataDir("G4INCLDATA")};
      G4String dataPathppbar(dataPath0 + "/rawppbarFS.dat");
      G4String dataPathnpbar(dataPath0 + "/rawnpbarFS.dat");
      G4String dataPathppbark(dataPath0 + "/rawppbarFSkaonic.dat");
      G4String dataPathnpbark(dataPath0 + "/rawnpbarFSkaonic.dat");
    #else
      Config const *theConfig=theNucleus->getStore()->getConfig();
      std::string path;
      if(theConfig)
        path = theConfig->getINCLXXDataFilePath();
      std::string dataPathppbar(path + "/rawppbarFS.dat");
      INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N ppbar final states" << dataPathppbar << '\n');
      std::string dataPathnpbar(path + "/rawnpbarFS.dat");
      INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N npbar final states" << dataPathnpbar << '\n');
      std::string dataPathppbark(path + "/rawppbarFSkaonic.dat");
      INCL_DEBUG("Reading https://doi.org/10.1016/j.physrep.2005.03.002 ppbar kaonic final states" << dataPathppbark << '\n');
      std::string dataPathnpbark(path + "/rawnpbarFSkaonic.dat");
      INCL_DEBUG("Reading https://doi.org/10.1007/BF02818764 and https://link.springer.com/article/10.1007/BF02754930 npbar kaonic final states" << dataPathnpbark << '\n');
   #endif
    /*std::string path = {"/home/zdemid/INCL/inclcode/data"};
    std::string dataPathppbar(path + "/rawppbarFS.dat");
    INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N ppbar final states" << dataPathppbar << '\n');
    std::string dataPathnpbar(path + "/rawnpbarFS.dat");
    INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N npbar final states" << dataPathnpbar << '\n');
    std::string dataPathppbark(path + "/rawppbarFSkaonic.dat");
    INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N ppbar kaonic final states" << dataPathppbark << '\n');
    std::string dataPathnpbark(path + "/rawnpbarFSkaonic.dat");
    INCL_DEBUG("Reading https://doi.org/10.1016/0375-9474(92)90362-N npbar kaonic final states" << dataPathnpbark << '\n');
    */
    //read probabilities and particle types from file
    std::vector<G4double> probabilities; //will store each FS yield
    std::vector<std::vector<std::string>> particle_types; //will store particle names
    G4double sum; //will contain a sum of probabilities of all FS in the file
    G4double kaonicFSprob=0.05; //probability to kave kaonic FS

    const G4bool isProton = ProtonIsTheVictim();
    G4int z = theNucleus->getZ(); //was modified in Cascade.cc
    G4int a = theNucleus->getA(); //was modified in Cascade.cc
    a++; //restoration of original A value before annihilation
    if(isProton == true){z++;} //restoration of original Z value before annihilation
    ThreeVector annihilationPosition;
    ParticleList starlist;
    ThreeVector mommy; //momentum to be assigned later

    //LETS GOOOOOOO!!!
    G4double rdm = Random::shoot(); 
    if(isProton == true){ //protonic annihilation
      INCL_DEBUG("Proton is the victim" << '\n');
      if(rdm < (1.-kaonicFSprob)){ // pionic FS was chosen
        INCL_DEBUG("pionic pp final state chosen" << '\n');
        sum = read_file(dataPathppbar, probabilities, particle_types);
        rdm = (rdm/(1.-kaonicFSprob))*sum; //99.88 normalize by the sum of probabilities in the file
        //now get the line number in the file where the FS particles are stored:
        G4int n = findStringNumber(rdm, probabilities)-1;
        if ( n < 0 ) return starlist;
        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
          if(particle_types[n][j] == "pi0"){
            Particle *p = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "pi-"){
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "pi+"){
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "omega"){
            Particle *p = new Particle(Omega, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "eta"){
            Particle *p = new Particle(Eta, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "rho-"){
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(pp);
          }
          else if(particle_types[n][j] == "rho+"){
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(pp);
          }
          else if(particle_types[n][j] == "rho0"){
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(pp);
          }
          else{
            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
            }
            std::cout << "Some non-existing FS particle detected when reading pbar FS files" << std::endl;
          }
        }
      }
      else{
        INCL_DEBUG("kaonic pp final state chosen" << '\n');
        sum = read_file(dataPathppbark, probabilities, particle_types);
        rdm = ((1-rdm)/kaonicFSprob)*sum;//2670 normalize by the sum of probabilities in the file
        //now get the line number in the file where the FS particles are stored:
        G4int n = findStringNumber(rdm, probabilities)-1;
        if ( n < 0 ) return starlist;
        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
          if(particle_types[n][j] == "pi0"){
            Particle *p = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "pi-"){
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "pi+"){
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "omega"){
            Particle *p = new Particle(Omega, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "eta"){
            Particle *p = new Particle(Eta, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "K-"){
            Particle *p = new Particle(KMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "K+"){
            Particle *p = new Particle(KPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "K0"){
            Particle *p = new Particle(KZero, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "K0b"){
            Particle *p = new Particle(KZeroBar, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else{
            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
            }
            std::cout << "Some non-existing FS particle detected when reading pbar FS files" << std::endl;
          }
        }
      }
    }
    else{ //neutronic annihilation
      INCL_DEBUG("Neutron is the victim" << '\n');
      if(rdm < (1.-kaonicFSprob)){ // pionic/kaonic choice
        INCL_DEBUG("pionic np final state chosen" << '\n');
        sum = read_file(dataPathnpbar, probabilities, particle_types);
        rdm = (rdm/(1.-kaonicFSprob))*sum; //99.95 normalize by the sum of probabilities in the file
        //now get the line number in the file where the FS particles are stored:
        G4int n = findStringNumber(rdm, probabilities)-1;
        if ( n < 0 ) return starlist;
        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
          if(particle_types[n][j] == "pi0"){
            Particle *p = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "pi-"){
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "pi+"){
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "omega"){
            Particle *p = new Particle(Omega, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "eta"){
            Particle *p = new Particle(Eta, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "rho-"){
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(pp);
          }
          else if(particle_types[n][j] == "rho+"){
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(pp);
          }
          else if(particle_types[n][j] == "rho0"){
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
            Particle *pp = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(pp);
          }
          else{
            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
            }
            std::cout << "Some non-existing FS particle detected when reading pbar FS files" << std::endl;
          }
        }
      }
      else{
        INCL_DEBUG("kaonic np final state chosen" << '\n');
        sum = read_file(dataPathnpbark, probabilities, particle_types);
        rdm = ((1-rdm)/kaonicFSprob)*sum;//3837 normalize by the sum of probabilities in the file
        //now get the line number in the file where the FS particles are stored:
        G4int n = findStringNumber(rdm, probabilities)-1;
        if ( n < 0 ) return starlist;
        for(G4int j = 0; j < static_cast<G4int>(particle_types[n].size()); j++){
          if(particle_types[n][j] == "pi0"){
            Particle *p = new Particle(PiZero, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "pi-"){
            Particle *p = new Particle(PiMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "pi+"){
            Particle *p = new Particle(PiPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "omega"){
            Particle *p = new Particle(Omega, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "eta"){
            Particle *p = new Particle(Eta, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "K-"){
            Particle *p = new Particle(KMinus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "K+"){
            Particle *p = new Particle(KPlus, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "K0"){
            Particle *p = new Particle(KZero, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else if(particle_types[n][j] == "K0b"){
            Particle *p = new Particle(KZeroBar, mommy, annihilationPosition);
            starlist.push_back(p);
          }
          else{
            INCL_ERROR("Some non-existing FS particle detected when reading pbar FS files");
            for(G4int jj = 0; jj < static_cast<G4int>(particle_types[n].size()); jj++){
              std::cout << "gotcha! " << particle_types[n][jj] << std::endl;
            }
            std::cout << "Some non-existing FS particle detected when reading pbar FS files" << std::endl;
          }
        }
      }
    }

    // Correction to the Q-value of the entering particle
    G4double theCorrection1 = theParticle->getEmissionPbarQvalueCorrection(a, z, isProton);
    G4double theCorrection2 = theParticle->getEmissionPbarQvalueCorrection(a, z, !isProton);
    G4double energyOfMesonStar;
    if(isProton == true){
      energyOfMesonStar = theParticle->getTableMass() + ParticleTable::getTableMass(a,z,0)
      -ParticleTable::getTableMass(a-1,z-1,0);
    }
    else{
      energyOfMesonStar = theParticle->getTableMass() + ParticleTable::getTableMass(a,z,0)
      -ParticleTable::getTableMass(a-1,z,0) + theCorrection2 - theCorrection1;
    }
    
    //compute energies of mesons with a phase-space model
    if(starlist.size() < 2){
      INCL_ERROR("should never happen, at least 2 final state particles!" << '\n');
    }
    else if(starlist.size() == 2){
      ParticleIter first = starlist.begin(); 
      ParticleIter last = std::next(first, 1); //starlist.end() gives an error of segfault, idk why
      G4double m1 = (*first)->getMass();
      G4double m2 = (*last)->getMass();
      G4double s = energyOfMesonStar*energyOfMesonStar;  
      G4double mom1 = std::sqrt(s/4 - (std::pow(m1,2) + std::pow(m2,2))/2 - std::pow(m1,2)*std::pow(m2,2)/s + (std::pow(m1,4) + 2*std::pow(m1*m2,2) + std::pow(m2,4))/(4*s));
      ThreeVector momentello = Random::normVector(mom1); //like raffaello :)
      (*first)->setMomentum(momentello);
      (*first)->adjustEnergyFromMomentum();
      (*last)->setMomentum(-momentello);
      (*last)->adjustEnergyFromMomentum();
      //std::cout << (*first)->getEnergy() << std::endl;
    }
    else{
      PhaseSpaceGenerator::generate(energyOfMesonStar, starlist);
      //ParticleIter first = starlist.begin(); 
      //std::cout << (*first)->getEnergy() << std::endl;
      //ParticleIter last = std::next(first, 1);
      //std::cout << (*last)->getEnergy() << std::endl;
    }
    
    return starlist;
  }  
   
  G4bool PbarAtrestEntryChannel::ProtonIsTheVictim(){
    if(theNucleus->getAnnihilationType() == PType){
      INCL_DEBUG("isProton" << '\n');
      return true; //proton is annihilated
    }
    else if(theNucleus->getAnnihilationType() == NType){
      INCL_DEBUG("isNeutron" << '\n');
      return false; //neutron is annihilated
    }
    else{
      INCL_ERROR("should never happen, n or p is your only choice!" << '\n');
      G4double rdm3 = Random::shoot(); 
      if(rdm3 >= 0.){  
        // it is set here for test
        return false;
      }
      else{
        return true;
      }
    }
  }

    //compute energy lost due to binding with electron shell
  G4double PbarAtrestEntryChannel::PbarCoulombicCascadeEnergy(G4int A, G4int Z){
    G4double N_ann = n_annihilation(A, Z);
    return ParticleTable::getINCLMass(antiProton)*(A/(A+1.))*(Z*Z/(N_ann*N_ann*2.*137.*137.)); 
    
  }
    /* Coulombic Cascade Energy formula in Bohr approximation taken from: 
     Precision spectroscopy of light exotic atoms D. Gotta
     Progress in Particle and Nuclear Physics 52 (2004) 133–195
     This is a crude approximation*/

  ThreeVector PbarAtrestEntryChannel::getAnnihilationPosition(){
    const G4bool isProton = ProtonIsTheVictim();
    G4int z = theNucleus->getZ(); //was modified in Cascade.cc
    G4int a = theNucleus->getA(); //was modified in Cascade.cc
    G4double n_ann = n_annihilation(a, z);
    a++; //not before the n_ann!

    if(isProton == true){z++;}
    G4double Rpmax = ParticleTable::getMaximumNuclearRadius(Proton, a, z);
    G4double Rnmax = ParticleTable::getMaximumNuclearRadius(Neutron, a, z);
    G4double probabilitymax = 0.; //the max value of the probability distribution
    G4double probability = 0.0; 
    G4double radius;

    //now we compute the max value of the probability distribution...
    if(isProton == true){

      for(radius = 0.0; radius < Rpmax; radius = radius + 0.001){  
        probability = overlapP(radius, n_ann);
              //INCL_WARN("radius, densityP, overlapP: " << radius << " "  << densityP(radius) << " " << probability << '\n');
        if(probability > probabilitymax)
        probabilitymax = probability; //now it should be the max value of overlapP function
      }
    }
    else{ //neutron

      for(radius = 0.0; radius < Rnmax; radius = radius + 0.001){     
        probability = overlapN(radius, n_ann);
              //INCL_WARN("radius, densityN, overlapN: " << radius << " "  << densityN(radius) << " " << probability << '\n');
        if(probability > probabilitymax)
        probabilitymax = probability; //now it should be the max value of overlapP function
      }
    }

    //we know the limits! start rejection algorithm!        
    G4double x = 0., y = 0.0001, p_for_x = 0.;
    G4double distance = 0.;
    if(isProton == true){  
      while(y >= p_for_x){
        x = Random::shoot() * Rpmax; // create uniformly random r
        y = Random::shoot() * probabilitymax; // create uniformly random prob
        p_for_x = overlapP(x, n_ann); //probability call for comparison
        if(y <= p_for_x){ //first cut-off is introduced for computational volume reduction
          distance = x;
        }
      }    
    }
    else{
      while(y >= p_for_x){
        x = Random::shoot() * Rnmax; // create uniformly random r
        y = Random::shoot() * probabilitymax; // create uniformly random prob
        p_for_x = overlapN(x, n_ann); //probability call for comparison
        if(y <= p_for_x){ //first cut-off is introduced for computational volume reduction
          distance = x;   
        }
      }
    }

    //FINAL POSITION VECTOR
    ThreeVector annihilationPosition(0., 0., -distance); //3D sphere of distance radius


    return annihilationPosition;
  }


  G4double PbarAtrestEntryChannel::n_annihilation(G4int A, G4int Z){
    const G4bool isProton = ProtonIsTheVictim();
    G4int z = Z; 
    G4int a = A; 
    a++;
    if(isProton == true){
      z++;
    }
    INCL_DEBUG("the original Z value is " << z << '\n');
    INCL_DEBUG("the original A value is " << a << '\n');
    G4double n_ann; //annihilation principal quantum number(interpolation from data H.Poth)
    if(z <= 1.){
      n_ann = 1.;
    }
    else if(z <= 4.){
      n_ann = 2.;
    }
    else if(z <= 11.){
      n_ann = 3.;      
    }
    else if(z <= 20.){
      n_ann = 4.;      
    }
    else if(z <= 32.){
      n_ann = 5.;      
    }
    else if(z <= 46.){
      n_ann = 6.;      
    }
    else if(z <= 61.){
      n_ann = 7.;      
    }
    else if(z <= 74.){
      n_ann = 8.;      
    }
    else if(z <= 84.){
      n_ann = 9.;      
    }
    else{
      n_ann = 10.;     
    }
    INCL_DEBUG("The following Pbar will annihilate with n = " << n_ann << '\n');

    return n_ann;
  }

  IAvatarList PbarAtrestEntryChannel::bringMesonStar(ParticleList const &pL, Nucleus * const n) {
    ThreeVector ann_position = getAnnihilationPosition();
    IAvatarList theAvatarList;
    for(ParticleIter p = pL.begin(), e = pL.end(); p!=e; ++p){
      (*p)->setPosition(ann_position);
      theAvatarList.push_back(new ParticleEntryAvatar(0.0, n, *p, APAR));
    }
    return theAvatarList;
  }

  void PbarAtrestEntryChannel::fillFinalState(FinalState *fs) {
    //const G4bool isProton = ProtonIsTheVictim();
    //G4int z = theNucleus->getZ(); //was modified in Cascade.cc
    //G4int a = theNucleus->getA(); //was modified in Cascade.cc
    //a++; //restoration of original A value before annihilation
    //if(isProton == true){z++;} //restoration of original Z value before annihilation
    const G4double energyBefore = theParticle->getEnergy(); 
    fs->addEnteringParticle(theParticle); 
    INCL_DEBUG("Entering particle added " << '\n');    
    fs->setTotalEnergyBeforeInteraction(energyBefore);
  }

}



