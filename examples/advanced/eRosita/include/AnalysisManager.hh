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
#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include "globals.hh"
#include <fstream>



class AnalysisManager {

 public:
   // The analysis class is designed to be a singleton (i.e. only one instance
   // can exist). A member function called Instance is defined, which allows 
   // the user to get a pointer to the existing instance or to create it if 
   // it does not yet exist
   // 
   static AnalysisManager* Instance();

   // The analysis class instance can be deleted by calling the Destroy
   // method (NOTE: The class destructor is protected, and can thus not be
   // called directly):
   static void Destroy(); 

//    // Member function used to score the total energy deposit
//   void ScoreTot(G4double eTot);

   // Member function used to dump hits
  void Score(G4double eDep);

 protected:
   // Constructor (protected):
   AnalysisManager();    

   // Destructor (protected): 
   virtual ~AnalysisManager();

   // Prevent copying
   AnalysisManager(const AnalysisManager& only);
   const AnalysisManager& operator=(const AnalysisManager& only);

private:
  // The static instance of the AnalysisManager class:
  static AnalysisManager* instance;
  std::ofstream outFile;
  std::ofstream outFileT;
};


#endif // ANALYSISMANAGER_HH
