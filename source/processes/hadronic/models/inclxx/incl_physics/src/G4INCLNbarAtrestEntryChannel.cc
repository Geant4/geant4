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

/*
 * G4INCLNbarAtrestEntryChannel.cc
 *
 *  \date Aug 9, 2024
 * \author Olivier Lourgo
 */
#ifdef INCLXX_IN_GEANT4_MODE
   #include "G4EnvironmentUtils.hh"
#endif
#include "G4INCLNbarAtrestEntryChannel.hh"
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

namespace G4INCL{

   NbarAtrestEntryChannel::NbarAtrestEntryChannel(Nucleus *n, Particle *p)
      :theNucleus(n), theParticle(p)
   {}
   NbarAtrestEntryChannel::~NbarAtrestEntryChannel(){}
   
   G4double NbarAtrestEntryChannel::read_file(std::string filename, std::vector<G4double>& probabilities, std::vector<std::vector<std::string>>& particle_types){
      std::ifstream file(filename);
      G4double sum_probs =0.0;
      if (file.is_open()){
         std::string line;
         while(getline(file,line)){
            std::istringstream iss(line);
            G4double prob;
            iss >> prob;
            sum_probs += prob;
            probabilities.push_back(prob);
            std::vector<std::string> types;
            std::string type;
            while (iss >>type){
               types.push_back(type);
            }
            particle_types.push_back(types);
         }
      }
      else std::cout << "ERROR no fread_file " << filename << std::endl;

      return sum_probs;
   }


   G4int NbarAtrestEntryChannel::findStringNumber(G4double rdm, std::vector<G4double> yields){
      G4int stringNumber =-1;
      G4double smallestsum =0.0;
      G4double biggestsum = yields[0];
      for (G4int i=0; i < static_cast<G4int>(yields.size() -1);i++){
         if (rdm >= smallestsum && rdm <= biggestsum){
            stringNumber = i+1;

         }
         smallestsum += yields[i];
         biggestsum += yields[i+1];
      }
      if (stringNumber==-1) stringNumber = static_cast<G4int>(yields.size());
      if (stringNumber==-1){
         INCL_ERROR("ERROR in findStringNumber (stringNumber=-1)");
         std::cout << "ERROR in findStringNumber" << std::endl;
      }
      return stringNumber;
   }

   G4double NbarAtrestEntryChannel::Pabs(G4double x, G4double value){
       const G4double r = value; // center of the gaussian
       const G4double sigma = 1;
       return std::exp(-std::pow(x-r,2)/(2*sigma*sigma));
    }

   G4double NbarAtrestEntryChannel::densityP(){ // return the r at which the gaussian of the interaction Probability(Pabs) is centered
      const G4bool isProton = ProtonIsTheVictim();
      G4int Z = theNucleus->getZ(); //was modified in Cascade.cc
      G4int A = theNucleus->getA(); //was modified in Cascade.cc
      G4double threshold_density = 0.10; //the maximum of the interaction probability is taken at 10% of maximum density 
      //https://doi.org/10.1016/0375-9474(82)90352-9 , Nuclear absorption of stopped antiprotons: Multipion-nucleus interactions, Iljinov, Nazaruk, Chigrinov
      A++; //restoration of original A value before annihilation
      if(isProton == true){Z++;} //restoration of original Z value before annihilation
      if(theNucleus->getAnnihilationType()==DNbarPPbarNType || theNucleus->getAnnihilationType()==DNbarNPbarPType || 
         theNucleus->getAnnihilationType()==DNbarNPbarNType || theNucleus->getAnnihilationType()==DNbarNPbarPType){A++;}
      if(theNucleus->getAnnihilationType()==DNbarPPbarPType || theNucleus->getAnnihilationType()==DNbarNPbarPType){Z++;}

      if(A > 19) {
         G4double radius = ParticleTable::getRadiusParameter(Proton, A, Z);
         G4double diffuseness = ParticleTable::getSurfaceDiffuseness(Proton, A, Z);
         G4double r_10 = diffuseness*std::log((1/threshold_density)-1) + radius;  //Radius for a Wood-Saxon
         return r_10;
      }else if(A <= 19 && A > 6) {
         G4double radius = ParticleTable::getRadiusParameter(Proton, A, Z);
         G4double diffuseness = ParticleTable::getSurfaceDiffuseness(Proton, A, Z);
         G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Proton, A, Z);
         NuclearDensityFunctions::ModifiedHarmonicOscillator rDensityFunction(radius, maximumRadius, diffuseness);
         //double r_10 = 0.01;
         //while (rDensityFunction(r_10)/(r_10*r_10) > threshold_density*rDensityFunction(0.01)/(0.01*0.01)) {
         // r_10 = r_10 + maximumRadius/100. ;
         //}
         G4double r_min = 0.01;
         G4double r_max = maximumRadius;
         G4double r_10 = (r_min + r_max)/2.;           
         while ((rDensityFunction(r_10)/(r_10*r_10) > 0.11*rDensityFunction(0.01)/(0.01*0.01)) || 
                (rDensityFunction(r_10)/(r_10*r_10) < 0.09*rDensityFunction(0.01)/(0.01*0.01))) {
           if (rDensityFunction(r_10)/(r_10*r_10) > 0.11*rDensityFunction(0.01)/(0.01*0.01)) {
             r_min = r_10;
           }
           else {
            r_max = r_10;
           }
           r_10 = (r_min + r_max)/2.;
         }
         return r_10;
      }else if(A <= 6 && A > 2) { // Gaussian distribution for light nuclei
         G4double radius = ParticleTable::getRadiusParameter(Proton, A, Z);
         G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Proton, A, Z);
         NuclearDensityFunctions::Gaussian rDensityFunction(maximumRadius, Math::oneOverSqrtThree * radius);
         //double r_10=std::sqrt(std::pow(Math::oneOverSqrtThree * radius,2)*std::log(2)); //start when the density is half the maximum
         //while (rDensityFunction(r_10)/(r_10*r_10) > threshold_density*rDensityFunction(0.01)/(0.01*0.01)) {
         // r_10 = r_10 + maximumRadius/500. ;
         //}
         G4double r_min = 0.01;
         G4double r_max = maximumRadius;
         G4double r_10 = (r_min + r_max)/2.;           
         while ((rDensityFunction(r_10)/(r_10*r_10) > 0.11*rDensityFunction(0.01)/(0.01*0.01)) || 
                (rDensityFunction(r_10)/(r_10*r_10) < 0.09*rDensityFunction(0.01)/(0.01*0.01))) {
           if (rDensityFunction(r_10)/(r_10*r_10) > 0.11*rDensityFunction(0.01)/(0.01*0.01)) {
             r_min = r_10;
           }
           else {
            r_max = r_10;
           }
           r_10 = (r_min + r_max)/2.;
         }
         return r_10;
      }else {
         INCL_ERROR("No nuclear density function for target A = "
               << A << " Z = " << Z << '\n');
         return 0.0;
      }

   }

   G4double NbarAtrestEntryChannel::densityN(){
      const G4bool isProton = ProtonIsTheVictim();
      G4int Z = theNucleus->getZ(); //was modified in Cascade.cc
      G4int A = theNucleus->getA(); //was modified in Cascade.cc
      G4double threshold_density = 0.10; //the maximum of the interaction probability is taken at 10% of maximum density
      //https://doi.org/10.1016/0375-9474(82)90352-9 , Nuclear absorption of stopped antiprotons: Multipion-nucleus interactions, Iljinov, Nazaruk, Chigrinov
      A++; //restoration of original A value before annihilation
      if(isProton == true){Z++;} //restoration of original Z value before annihilation
      if(theNucleus->getAnnihilationType()==DNbarPPbarNType || theNucleus->getAnnihilationType()==DNbarNPbarPType || 
         theNucleus->getAnnihilationType()==DNbarNPbarNType || theNucleus->getAnnihilationType()==DNbarNPbarPType){A++;}
      if(theNucleus->getAnnihilationType()==DNbarPPbarPType || theNucleus->getAnnihilationType()==DNbarNPbarPType){Z++;}

      if(A > 19) {
         G4double radius = ParticleTable::getRadiusParameter(Neutron, A, Z);
         G4double diffuseness = ParticleTable::getSurfaceDiffuseness(Neutron, A, Z);
         G4double r_10 = diffuseness*std::log((1/threshold_density)-1) + radius;  //Radius for a Wood-Saxon
         return r_10;
      } else if(A <= 19 && A > 6) {
         G4double radius = ParticleTable::getRadiusParameter(Neutron, A, Z);
         G4double diffuseness = ParticleTable::getSurfaceDiffuseness(Neutron, A, Z);
         G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Neutron, A, Z);
         NuclearDensityFunctions::ModifiedHarmonicOscillator rDensityFunction(radius, maximumRadius, diffuseness);
         G4double r_10 = 0.01;
         while (rDensityFunction(r_10)/(r_10*r_10) > threshold_density*rDensityFunction(0.01)/(0.01*0.01)) {
          r_10 = r_10 + maximumRadius/100. ;
         }
         return r_10;
      } else if(A <= 6 && A > 2) { // Gaussian distribution for light nuclei
         G4double radius = ParticleTable::getRadiusParameter(Neutron, A, Z);
         G4double maximumRadius = ParticleTable::getMaximumNuclearRadius(Neutron, A, Z);
         NuclearDensityFunctions::Gaussian rDensityFunction(maximumRadius, Math::oneOverSqrtThree * radius);
         G4double r_10=std::sqrt(std::pow(Math::oneOverSqrtThree * radius,2)*std::log(2)); //start when the density is half the maximum
         while (rDensityFunction(r_10)/(r_10*r_10) > threshold_density*rDensityFunction(0.01)/(0.01*0.01)) {
          r_10 = r_10 + maximumRadius/500. ;
         }
         return r_10;
      } else {
         INCL_ERROR("No nuclear density function for target A = "
               << A << " Z = " << Z << '\n');
         return 0.0;
      }


   }

   G4double NbarAtrestEntryChannel::overlapN(G4double &x){
      return Pabs(x,densityN());
   }
   G4double NbarAtrestEntryChannel::overlapP(G4double &x){
      return Pabs(x,densityP());

   }
   ParticleList NbarAtrestEntryChannel::makeMesonStar() {//This function creates a set of mesons with momenta
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
      G4String dataPathnbarp(dataPath0 + "/rawnbarpFS.dat");
      G4String dataPathnbarn(dataPath0 + "/rawnbarnFS.dat");
      G4String dataPathnbarnk(dataPath0 + "/rawppbarFSkaonic.dat");
      G4String dataPathnbarpk(dataPath0 + "/rawnbarpFSkaonic.dat");
    #else
      Config const *theConfig=theNucleus->getStore()->getConfig();
      std::string path;
      if(theConfig)
        path = theConfig->getINCLXXDataFilePath();
      std::string dataPathnbarn(path + "/rawnbarnFS.dat");
      INCL_DEBUG("Reading nbarn final states" << dataPathnbarn << '\n');
      std::string dataPathnbarp(path + "/rawnbarpFS.dat");
      INCL_DEBUG("Reading nbarp final states" << dataPathnbarp << '\n');
      std::string dataPathnbarnk(path + "/rawppbarFSkaonic.dat");
      INCL_DEBUG("Reading nbarn kaonic final states" << dataPathnbarnk << '\n');
      std::string dataPathnbarpk(path + "/rawnbarpFSkaonic.dat");
      INCL_DEBUG("Reading nbarp kaonic final states" << dataPathnbarpk << '\n');
   #endif
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
        sum = read_file(dataPathnbarp, probabilities, particle_types);
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
        sum = read_file(dataPathnbarpk, probabilities, particle_types);
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
        sum = read_file(dataPathnbarn, probabilities, particle_types);
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
        sum = read_file(dataPathnbarnk, probabilities, particle_types);
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
    G4int stra = theNucleus->getS(); 
    G4double energyOfMesonStar;
    if(theNucleus->isNucleusNucleusCollision()==false){//antiNeutron
    if(isProton == true){
      energyOfMesonStar = theParticle->getEnergy() + ParticleTable::getTableMass(a,z,stra)
                         -ParticleTable::getTableMass(a-1,z,stra);
    }
    else{
      energyOfMesonStar = theParticle->getEnergy() + ParticleTable::getTableMass(a,z,stra)
                         -ParticleTable::getTableMass(a-1,z,stra);
    }
    } else if(theNucleus->isNucleusNucleusCollision()==true){//antiComposite : job is done in the Antinuclei file
          return starlist;
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
   
   G4bool NbarAtrestEntryChannel::ProtonIsTheVictim(){
      if(theNucleus->getAnnihilationType() == PType || theNucleus->getAnnihilationType() == DNbarPPbarPType || theNucleus->getAnnihilationType() == DNbarPPbarNType ){
         return true; //a proton is annihilated
      }
      else if(theNucleus->getAnnihilationType() == NType || theNucleus->getAnnihilationType() == DNbarNPbarPType || theNucleus->getAnnihilationType() == DNbarNPbarNType){
         return false; // a neutron is annihilated
      }
      else{
         INCL_ERROR("should never happen, n or p is your only choise" << '\n');
         G4double rdm3 = Random::shoot();
         if(rdm3 >= 0.){
            return false;
         }
         else{
            return true;
         }
      }
   }
   ThreeVector NbarAtrestEntryChannel::getAnnihilationPosition(){
    const G4bool isProton = ProtonIsTheVictim();
    G4int z = theNucleus->getZ(); //was modified in Cascade.cc
    G4int a = theNucleus->getA(); //was modified in Cascade.cc
    a++; 

    if(isProton == true){z++;}
    if(theNucleus->getAnnihilationType()==DNbarPPbarNType || theNucleus->getAnnihilationType()==DNbarNPbarPType || 
       theNucleus->getAnnihilationType()==DNbarNPbarNType || theNucleus->getAnnihilationType()==DNbarNPbarPType){a++;}
    if(theNucleus->getAnnihilationType()==DNbarPPbarPType || theNucleus->getAnnihilationType()==DNbarNPbarPType){z++;}
    G4double Rpmax = ParticleTable::getMaximumNuclearRadius(Proton, a, z);
    G4double Rnmax = ParticleTable::getMaximumNuclearRadius(Neutron, a, z);
    G4double probabilitymax = 0.; //the max value of the probability distribution
    G4double probability = 0.0; 
    G4double radius;

    //now we compute the max value of the probability distribution...
    if(isProton == true){

      for(radius = 0.0; radius < Rpmax; radius = radius + 0.001){  
        probability = overlapP(radius);
              //INCL_WARN("radius, densityP, overlapP: " << radius << " "  << densityP(radius) << " " << probability << '\n');
        if(probability > probabilitymax)
        probabilitymax = probability; //now it should be the max value of overlapP function
      }
    }
    else{ //neutron

      for(radius = 0.0; radius < Rnmax; radius = radius + 0.001){     
        probability = overlapN(radius);
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
        p_for_x = overlapP(x); //probability call for comparison
        if(y <= p_for_x){ //first cut-off is introduced for computational volume reduction
          distance = x;
        }
      }    
    }
    else{
      while(y >= p_for_x){
        x = Random::shoot() * Rnmax; // create uniformly random r
        y = Random::shoot() * probabilitymax; // create uniformly random prob
        p_for_x = overlapN(x); //probability call for comparison
        if(y <= p_for_x){ //first cut-off is introduced for computational volume reduction
          distance = x;   
        }
      }
    }

    //FINAL POSITION VECTOR
    //ThreeVector annihilationPosition(0., 0., -distance); //3D sphere of distance radius
    G4double ctheta = (1.-2.*Random::shoot());
    G4double stheta = std::sqrt(1.-ctheta*ctheta);
    G4double phi = Math::twoPi*Random::shoot();
    ThreeVector annihilationPosition(distance*stheta * std::cos(phi), distance*stheta * std::sin(phi), distance*ctheta); //3D sphere of distance radius

    return annihilationPosition;
  }

   IAvatarList NbarAtrestEntryChannel::bringMesonStar(ParticleList const &pL, Nucleus * const n) {
    ThreeVector ann_position = getAnnihilationPosition();
    IAvatarList theAvatarList;
    for(ParticleIter p = pL.begin(), e = pL.end(); p!=e; ++p){
      (*p)->setPosition(ann_position);
      theAvatarList.push_back(new ParticleEntryAvatar(0.0, n, *p, ANAR));
    }
    return theAvatarList;
  }
   void NbarAtrestEntryChannel::fillFinalState(FinalState *fs) {
    //const bool isProton = ProtonIsTheVictim();
    //int z = theNucleus->getZ(); //was modified in Cascade.cc
    //int a = theNucleus->getA(); //was modified in Cascade.cc
    //a++; //restoration of original A value before annihilation
    //if(isProton == true){z++;} //restoration of original Z value before annihilation
    const G4double energyBefore = theParticle->getEnergy(); 
    fs->addEnteringParticle(theParticle); 
    INCL_DEBUG("Entering particle added " << '\n');    
    fs->setTotalEnergyBeforeInteraction(energyBefore);
  }

}