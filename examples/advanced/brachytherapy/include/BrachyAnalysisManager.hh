
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

class ITree;
class IHistogramFactory;
class IAnalysisFactory;
class ITupleFactory;
class ITuple;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class BrachyAnalysisManager
{
public:
  
   ~BrachyAnalysisManager();

  void book();
  
  void finish();

  static BrachyAnalysisManager* getInstance();
  
  void analyse(G4double,G4double,G4float);
  void  hist(G4double,G4double,G4float);
  void Spectrum(G4double);

private:

  G4double xx,zz;
  G4float  en; 
  G4double  x,z;
 static BrachyAnalysisManager* instance;

private:
   BrachyAnalysisManager();

private:

  IAnalysisFactory  *aFact;
  ITree             *theTree;
  IHistogramFactory *histFact;
  ITupleFactory     *tupFact;

};

#endif




