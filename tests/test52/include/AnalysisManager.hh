#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include <string>
#include "DataManager.hh"


class AnalysisManager : public DataManager {

 public:
   static AnalysisManager* Instance(double zLow = 0.0, 
                                    double zUp = 10.0);
   static void Destroy(); 
 
   void CreateCalorimeter(double pos, 
                          double thickn, 
                          double rad);
   void ScoreEnergyDeposit(double en, 
                           double x, 
                           double y, 
                           double z,
                           std::string type="");
   void ScoreParticleEnergy(double en, 
                            double x, 
                            double y, 
                            double z,
                            std::string type);
   void ScoreEnteringParticles(double en, 
                               double x, 
                               double y, 
                               double z,
                               std::string type);
   void ScoreExitingParticles(double en, 
                              double x, 
                              double y, 
                              double z,
                              std::string type, 
                              int id);
   void PrintResults();

 protected: 
   AnalysisManager(double zLow, double zUp);
   ~AnalysisManager();

 private:
   static AnalysisManager* instance;
   AIDA::IHistogram1D* energyInCalorimeters;
   std::vector<double> binBoundaries;

   double totEnergyDeposit;
   double totEnergyDepositSumSquares;

   double elecEnergyEnterTarget;
   double photEnergyExitTarget;
   double photEnergyExitTargetSumSquares;
   double elecEnergyExitTarget;
   double elecEnergyExitTargetSumSquares;
   double primElecEnergyExitTarget;

   int nmbEnergyDeposits;
   int nmbElecEnterTarget;
   int nmbElecExitTarget;
   int nmbPrimElecExitTarget;
   int nmbPhotExitTarget;
};

#endif // ANALYSISMANAGER_HH
