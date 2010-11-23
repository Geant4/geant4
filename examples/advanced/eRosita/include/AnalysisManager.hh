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

   // Member function used to score the total energy deposit
  void ScoreTot(G4double eTot);

   // Member function used to dump hits
  void Score(G4int id, G4double eDep, G4double x, G4double y, G4double z);

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
