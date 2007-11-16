#ifndef CALDATACOLLECTOR_HH
#define CALDATACOLLECTOR_HH

#include <string>
#include "DataManager.hh"


class CalDataCollector : public DataManager {

 public:
   CalDataCollector(double zLow,double zUp,double rad);
   ~CalDataCollector();

   void ScoreEnergyDeposit(double en, 
                           double x, 
                           double y, 
                           double z,
                           std::string type="");
   void ScoreParticleEnergy(double en, 
                            double x, 
                            double y, 
                            double z,
                            std::string type) {}
   void PrintResults();

 private:
   double radius;

   double totEnergyDeposit;
   double totEnergyDepositSumSquares;
   double elecEnergyDeposit;
   double photEnergyDeposit;
   int    nmbEnergyDeposits;

   AIDA::IHistogram2D* latEnergyDeposit;
};

#endif // CALDATACOLLECTOR_HH
