
#include <iostream>
#include <fstream>
#include <cmath>
#include <AIDA/AIDA.h>
#include "AnalysisManager.hh"
#include "CalDataCollector.hh"


AnalysisManager* AnalysisManager::instance = 0;


AnalysisManager* AnalysisManager::Instance(double zLow, double zUp) {
 
  if(instance == 0) {
    instance = new AnalysisManager(zLow, zUp);
  }
  
  return instance;
}


void AnalysisManager::Destroy() {

  if(!instance == 0) {

     delete instance;
     instance = 0;
  }
}


AnalysisManager::AnalysisManager(double zLow, double zUp) :
  DataManager(zLow,zUp),
  energyInCalorimeters(0),
  totEnergyDeposit(0),
  totEnergyDepositSumSquares(0),
  elecEnergyEnterTarget(0),
  photEnergyExitTarget(0),
  photEnergyExitTargetSumSquares(0),
  elecEnergyExitTarget(0),
  elecEnergyExitTargetSumSquares(0),
  primElecEnergyExitTarget(0),
  nmbEnergyDeposits(0),
  nmbElecEnterTarget(0),
  nmbElecExitTarget(0),
  nmbPrimElecExitTarget(0),
  nmbPhotExitTarget(0) {

}


AnalysisManager::~AnalysisManager() {

}


void AnalysisManager::CreateCalorimeter(double pos, 
                                        double thickn, 
                                        double rad) {

  double zLow = pos - 0.5 * thickn;
  double zUp  = pos + 0.5 * thickn;
  
  AddDataCollector(new CalDataCollector(zLow,zUp,rad));

  binBoundaries.push_back(zLow);
  sort(binBoundaries.begin(),binBoundaries.end());
}


void AnalysisManager::ScoreParticleEnergy(double en, 
                                          double x, 
                                          double y, 
                                          double z,
                                          std::string ptype) {

  if(z >= GetLowerBound() && z < GetUpperBound()) {
     DataManager::ScoreParticleEnergy(en,x,y,z,ptype);
  }
}


void AnalysisManager::ScoreEnergyDeposit(double en, 
                                         double x, 
                                         double y, 
	                                 double z,
                                         std::string ptype) {
 
  totEnergyDeposit += en;
  totEnergyDepositSumSquares += en * en;
  nmbEnergyDeposits++;
  DataManager::ScoreEnergyDeposit(en,x,y,z,ptype);

  if(energyInCalorimeters == 0 && binBoundaries.size() > 1) {

     std::string histName = "EnergyVsDepth";
     std::string histTitle = "Energy deposit in calorimeter slices"; 

     binBoundaries.push_back(GetUpperBound());

     energyInCalorimeters =
         histogramFactory -> createHistogram1D(histName,
                                               histTitle,
                                               binBoundaries); 
  }

  if(energyInCalorimeters != 0) energyInCalorimeters -> fill(z, en);
}


void AnalysisManager::ScoreEnteringParticles(double en, 
                                             double x, 
                                             double y, 
                                             double z,
                                             std::string ptype) {

  if(ptype == "electron" || ptype == "e-") {
     elecEnergyEnterTarget += en;
     nmbElecEnterTarget++;
  }
}


void AnalysisManager::ScoreExitingParticles(double en, 
                                            double x, 
                                            double y, 
	                                    double z,
                                            std::string ptype,
                                            int id) {
 
  if(ptype == "electron" || ptype == "e-") {
     elecEnergyExitTarget += en;     
          
     if(en >= 0.0) { 
        nmbElecExitTarget += 1;
        elecEnergyExitTargetSumSquares += en * en;
     }
     if(en < 0.0) {
        nmbElecExitTarget -= 1;
        elecEnergyExitTargetSumSquares -= en * en;
     }

     if(id == 1) {
        primElecEnergyExitTarget += en;

        if(en >= 0.0) nmbPrimElecExitTarget += 1;
        if(en < 0.0) nmbPrimElecExitTarget -= 1;
     }
  }
  if(ptype == "gamma") {
     photEnergyExitTarget += en;     
     photEnergyExitTargetSumSquares += en * en;     
     nmbPhotExitTarget += 1;
  }
}


void AnalysisManager::PrintResults() {

  os() << "--------------------------------------------------------" 
       << std::endl;
  os() << "TARGET: " 
       << GetLowerBound() << " to "
       << GetUpperBound() << std::endl;
  os() << "--------------------------------------------------------" 
       << std::endl;

  double totEnergyDepositErrSqu =         
           totEnergyDepositSumSquares / (totEnergyDeposit * 
           totEnergyDeposit) - 1.0 / double(nmbEnergyDeposits);  
  double totEnergyDepositRelErr = sqrt(totEnergyDepositErrSqu);        

  double elecEnergyExitTargetErrSqu =
           elecEnergyExitTargetSumSquares / (elecEnergyExitTarget *
	   elecEnergyExitTarget) - 1.0 / double(nmbElecExitTarget);
  double elecEnergyExitTargetRelErr = sqrt(elecEnergyExitTargetErrSqu);
 
  double photEnergyExitTargetErrSqu =
           photEnergyExitTargetSumSquares / (photEnergyExitTarget *
	   photEnergyExitTarget) - 1.0 / double(nmbPhotExitTarget);
  double photEnergyExitTargetRelErr = sqrt(photEnergyExitTargetErrSqu);


  os() << "  Total Energy Deposit: " 
       << totEnergyDeposit 
       << "   Rel. Error: "
       << totEnergyDepositRelErr
       << std::endl
       << "  Electron Energy/Number Entering Target: " 
       << elecEnergyEnterTarget << " "
       << nmbElecEnterTarget
       << std::endl            
       << "  Electron Energy/Number Exiting Target:  " 
       << elecEnergyExitTarget << " "
       << nmbElecExitTarget 
       << "  Rel. Error Energy: "
       << elecEnergyExitTargetRelErr
       << std::endl            
       << "  Primary Electron Energy/Number Exiting Target:  " 
       << primElecEnergyExitTarget << " "
       << nmbPrimElecExitTarget
       << std::endl            
       << "  Photon Energy/Number Exiting Target:    " 
       << photEnergyExitTarget << " "
       << nmbPhotExitTarget
       << "  Rel. Error Energy: "
       << photEnergyExitTargetRelErr
       << std::endl;            
 
  double backScEnergy = 
       1.0 - ((totEnergyDeposit+photEnergyExitTarget)/elecEnergyEnterTarget); 

  os() << "  Fraction of incid. electr. energy escaping as photons: "
       << photEnergyExitTarget/elecEnergyEnterTarget 
       << std::endl
       << "  Fraction of electron energy backscattered (indir. calcul.): "
       << backScEnergy
       << std::endl
       << "  Fraction of electron energy backscattered (direct calcul.): "
       << elecEnergyExitTarget/elecEnergyEnterTarget
       << std::endl
       << "  Fraction of electron energy backscattered, only primaries (direct calcul.): "
       << primElecEnergyExitTarget/elecEnergyEnterTarget
       << std::endl
       << "  Fraction of electrons backscattered (direct calcul.): "
       << double(nmbElecExitTarget)/double(nmbElecEnterTarget)
       << std::endl     
       << "  Fraction of electrons backscattered, only primaries " 
       << "(direct calcul.): "
       << double(nmbPrimElecExitTarget)/double(nmbElecEnterTarget)
       << std::endl;      
       
  DataManager::PrintResults();
  os() << "--------------------------------------------------------";

}      

