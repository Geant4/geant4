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
#ifndef DATAMANAGER_HH
#define DATAMANAGER_HH

#include <string>
#include <vector>
#include <fstream>
#include <utility>

class StorageManager;
class DataManager;
namespace AIDA {
   class IHistogramFactory;
   class IHistogram1D;
   class IHistogram2D;
}

typedef std::pair<double,DataManager*> coll;
typedef std::vector<coll> collector;


class zCompare {

 public:
   bool operator() (const coll& l,const coll& r) const 
                              { return keyLess(l.first,r.first); }

   bool operator() (const coll& l,const coll::first_type& k) const 
                              { return keyLess(l.first,k); }

   bool operator() (const coll::first_type& k,const coll& r) const 
                              { return keyLess(k,r.first); }

 private:
    bool keyLess(const coll::first_type& k1,
                 const coll::first_type& k2) const
	                       { return k1 < k2; }
};


class DataManager {
 
 public:
   DataManager(double zLow, double zUp);
   virtual ~DataManager();
   
   virtual void ScoreEnergyDeposit(double en, 
                                   double x, 
                                   double y, 
                                   double z,
                                   std::string type="all");
   virtual void ScoreParticleEnergy(double en, 
                                    double x, 
                                    double y, 
                                    double z,
                                    std::string type);
   virtual void PrintResults();

   void AddDataCollector(DataManager* comp); 
   bool HasDataCollector(double z);
   bool OverlapsWithOtherChild(double zLow, double zUp);
   bool IsContained(double zLow, double zUp);
   DataManager* MatchingChildDataCollector(double zLow, double zUp);

   double GetUpperBound() { return zAxisUpperBound; }
   double GetLowerBound() { return zAxisLowerBound; }

 protected:
   AIDA::IHistogramFactory* histogramFactory;

   std::ostream& os();

 private:
   StorageManager* storageManager;
   collector dataCollectors;
   collector garbage;

   double zAxisLowerBound;
   double zAxisUpperBound;   
   double eps;
};

#endif // DATAMANAGER_HH
