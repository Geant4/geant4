
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IAnalysisFactory.h"
namespace AIDA{
class ITree;
class IHistogramFactory;
class IAnalysisFactory;
class ITupleFactory;
class ITuple;
class ITreeFactory;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class BrachyAnalysisManager
{
public:

~BrachyAnalysisManager();

void book();

void finish();

static BrachyAnalysisManager* getInstance();

void fill_Tuple(G4double,G4double,G4double,G4float);
void hist(G4double,G4double,G4float);
void Spectrum(G4double);

private:

G4double xx,zz,yy;
G4float  en; 
G4double  x,y,z;
static BrachyAnalysisManager* instance;

private:
BrachyAnalysisManager();

private:

AIDA::IAnalysisFactory*  aFact;
AIDA::ITree*             theTree;
AIDA::IHistogramFactory *histFact;
AIDA::ITupleFactory     *tupFact;
AIDA::ITreeFactory      *treeFact;
AIDA::IHistogram2D *h1;
AIDA::IHistogram1D *h2;
AIDA::ITuple *ntuple;
};

#endif




