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
#ifndef STORAGEMANAGER_HH
#define STORAGEMANAGER_HH

#include <string>
#include <fstream>

namespace AIDA {
   class IAnalysisFactory;
   class ITree;
   class IHistogramFactory;
   class IHistogram1D;
   class IHistogram2D;
}


class StorageManager {

 public:
   static StorageManager* Instance(std::string fileNameBase);
   static void Destroy(); 

   AIDA::IHistogramFactory* GetIHistogramFactory();
   std::ostream& GetOutputFileStream();

 protected:
   virtual ~StorageManager();
   StorageManager(std::string fileNameBase);

   StorageManager(const StorageManager& only);
   const StorageManager& operator=(const StorageManager& only);

 private:
   static StorageManager* instance;

   AIDA::IHistogramFactory* histogramFactory;
   AIDA::IAnalysisFactory* analysisFactory;
   AIDA::ITree* tree;

   std::ofstream os;
};

#endif // STORAGEMANAGER_HH
