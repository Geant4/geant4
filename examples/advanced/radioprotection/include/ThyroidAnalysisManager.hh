
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

class ThyroidAnalysisManager
{
public:
 

   ~ThyroidAnalysisManager();

  void book();
  
  void finish();

  static ThyroidAnalysisManager* getInstance();
  
  void analyse(G4double,G4double,G4double,G4float);
  void  hist(G4double,G4double,G4float);
  void Spectrum(G4double);

private:

  G4double xx,yy,zz;
  G4float  en; 
  G4double  x,y,z;
 static ThyroidAnalysisManager* instance;

private:
   ThyroidAnalysisManager();

private:

  IAnalysisFactory  *aFact;
  ITree             *theTree;
  IHistogramFactory *histFact;
  ITupleFactory     *tupFact;

};

#endif




