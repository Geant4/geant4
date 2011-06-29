
#include <sstream>
#include <fstream>
#include <cmath>
#include <AIDA/AIDA.h>
#include "CalDataCollector.hh"


CalDataCollector::CalDataCollector(double zLow,double zUp,double rad) :
    DataManager(zLow,zUp),
    radius(rad),
    totEnergyDeposit(0),
    totEnergyDepositSumSquares(0),
    elecEnergyDeposit(0),
    photEnergyDeposit(0),
    nmbEnergyDeposits(0),
    latEnergyDeposit(0) {

}


CalDataCollector::~CalDataCollector() {

}


void CalDataCollector::ScoreEnergyDeposit(double en, 
                                          double x, 
                                          double y, 
					  double z,
                                          std::string ptype) {

  if(latEnergyDeposit == 0) {

     double center = 
              GetLowerBound() + (GetUpperBound() - GetLowerBound()) * 0.5;
     std::stringstream s;
     s << center;
     std::string histName = "z=" + s.str();
     std::string histTitle = "Lateral energy distr. in calorimeter at " 
                          + histName;

     latEnergyDeposit =
         histogramFactory -> createHistogram2D(histName,
                                               histTitle,
                                               50,-radius,radius,
                                               50,-radius,radius); 
  }

  if(z >= GetLowerBound() && z < GetUpperBound()) {

     if(ptype == "electron" || ptype == "e-")  elecEnergyDeposit += en;
     if(ptype == "gamma")     photEnergyDeposit += en;

     totEnergyDeposit += en;
     totEnergyDepositSumSquares += en * en;
     nmbEnergyDeposits++;
     latEnergyDeposit -> fill(x, y, en);
     
     DataManager::ScoreEnergyDeposit(en, x, y, z, ptype);
  }
}


void CalDataCollector::PrintResults() {

  os() << "--------------------------------------------------------" 
       << std::endl;
  os() << "CALORIMETER: " 
       << GetLowerBound() << " to "
       << GetUpperBound() << std::endl;
  os() << "--------------------------------------------------------" 
       << std::endl;
 
  double totEnergyDepositErrSqu =         
           totEnergyDepositSumSquares / (totEnergyDeposit * 
           totEnergyDeposit) - 1.0 / double(nmbEnergyDeposits);  
  double totEnergyDepositRelErr = std::sqrt(totEnergyDepositErrSqu);  

  double thickness = GetUpperBound() - GetLowerBound();
  double center = GetLowerBound() + 0.5 * thickness;


  os() << "Center and Thickness:  "
       << center << " " 
       << thickness 
       << "   Total Energy Deposit: "
       << totEnergyDeposit 
       << "   Rel. Error: "
       << totEnergyDepositRelErr
       << std::endl
       << "  Electron Energy Deposit: " 
       << elecEnergyDeposit
       << std::endl            
       << "  Photon Energy Deposit:   " 
       << photEnergyDeposit
       << std::endl;            

  DataManager::PrintResults();
}
