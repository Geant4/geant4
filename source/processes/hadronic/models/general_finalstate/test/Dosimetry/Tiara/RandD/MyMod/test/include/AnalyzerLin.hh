#ifndef __SYSANALYZER_DEFINED__
#define __SYSANALYZER_DEFINED__

#include "globals.hh"

class IVector;

enum HistoOptions {henGRID=0,henDVY,henDVX,henSTAT,henTIC,henBOX,henLINX,henLINY,xenAST};

#define AA_NONE 0
#define AA_BOX  1
#define AA_TENT 2

void alSetAA(int AA);
G4int alAddHisto(G4String szName,G4int nBins,G4double dMin,G4double dMax);
void alDeleteHisto(G4int nHisto);
void alShowHisto(G4int nHisto,G4int Style=0);
void alHideHisto(G4int nHisto);
void alMakePS(G4int nHisto,G4String file);
void alOpenStorage(G4String szName);
void alSaveHisto(G4int nHisto);
void alRemoveFromStore(G4int nHisto);
void alInit();
void alAddData(G4int nHist,IVector* pVecor);
void alAdditionalHisto(G4int nHisto,const char* Storage,const char* label);
void alAddDataFromFile(G4int nHist,G4String szFileName);
void ShowViewer(G4bool bYes);
void alSetProperty(G4int nHisto,G4String szKey,G4String szValue);
void alReset();
void alClear();
void alFill(G4int nHisto,G4double dValue,G4double weight,G4double EvWeight=1,bool bCalcErrors=true);
void alScaleAllHistos(G4double val);
void alSubstract(G4int nHisto,G4int nHistoOut,G4double dMaxEnergy,G4double dMinEnergy=0.);
void alCreateSmooth(G4int nHisto,G4int nTargetHisto,G4double dMinEnergy);
double alGetMinimumEntries(int LastID);
unsigned alGetEntries(int nHisto);
float alCalculateIntegral(int nHisto);
void alAddEvent(int nHisto);
void alRefreshHisto(bool bAddZone);
void alSetTitle(char* title);
void alListHistos();
void alAddTitles(int nHisto,char* XTitle,char* YTitle);
void alSetBit(int nHisto,unsigned Option,bool bSet);
void alReadOptions(unsigned int nHist,char* FileName);
void alDeleteOptions(unsigned int nHist,G4String szKey);
void alListOptions(unsigned int nHisto);
bool alRebin(unsigned nHistoFirst,unsigned nHistoLast,unsigned bias=0);
double alGetEnergy(int nHisto,double Energy);
double alGetEnergyRange(int nHisto,double Energy);
#endif
