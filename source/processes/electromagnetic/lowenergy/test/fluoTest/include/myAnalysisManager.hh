//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifdef G4ANALYSIS_USE
#ifndef myAnalysisManager_h
#define myAnalysisManager_h 1

#include "G4VAnalysisManager.hh"

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

class myAnalysisMessenger;
class myDetectorConstruction;
class IHistogramFactory;
class IHistogram1D;
class IPlotter;
class IVectorFactory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myAnalysisManager: public G4VAnalysisManager
{
public:
  myAnalysisManager(myDetectorConstruction*);
  virtual ~myAnalysisManager();
  
public:
  G4bool RegisterAnalysisSystem(G4VAnalysisSystem*);
  IHistogramFactory* GetHistogramFactory(const G4String&);
 
  void Store(IHistogram* = 0, const G4String& = "");
  void Plot(IHistogram* = 0);
  void InsertKEnergy(double gKe);
  void InsertIoniEnergy(double Ioe);
  void InsertPhotoEnergy(double Phe);
  void InsertBremEnergy(double Bre);
  void InsertComptEnergy(double Coe);
  void InsertConvEnergy(double Cve);
  void InsertRaylEnergy(double Rae);
  void InsertOutEnergy(double Goe);
 
  void InserteKEnergy(double eKe);
  void InserteIoniEnergy(double eIoe);
  void InsertePhotoEnergy(double ePhe);
  void InserteBremEnergy(double eBre);
  void InserteComptEnergy(double eCoe);
  void InserteConvEnergy(double eCve);
  void InserteRaylEnergy(double eRae);
  void InserteOutEnergy(double Eoe);

 void BeginOfRun();
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);

  void SetHisto1DDraw(G4String str) {histo1DDraw = str;};
  void SetHisto1DSave(G4String str) {histo1DSave = str;};
 
private:
  G4VAnalysisSystem* analysisSystem;
  IPlotter* pl;
  IVectorFactory* fVectorFactory;
  IHistogramFactory* histoFactory;

  IHistogram1D*  histoGammaKenergy;
  IHistogram1D* histoIoniEnergy;
  IHistogram1D*  histoPhotoEnergy;
  IHistogram1D* histoBremEnergy;
  IHistogram1D*  histoComptEnergy;
  IHistogram1D*  histoConvEnergy;
  IHistogram1D*  histoRaylEnergy;
  IHistogram1D*  histoGammaOutEnergy;
  IHistogram1D*  histoElecKenergy;
  IHistogram1D* histoeIoniEnergy;
  IHistogram1D*  histoePhotoEnergy;
  IHistogram1D* histoeBremEnergy;
  IHistogram1D*  histoeComptEnergy;
  IHistogram1D*  histoeConvEnergy;
  IHistogram1D*  histoeRaylEnergy;
  IHistogram1D*  histoElecOutEnergy;
  myDetectorConstruction*    Detector;

  G4String histo1DDraw;
  G4String histo1DSave;
 

  myAnalysisMessenger* analysisMessenger;
};


#endif
#endif



