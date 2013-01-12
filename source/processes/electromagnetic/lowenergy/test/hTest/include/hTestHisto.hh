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
#ifndef hTestHisto_h
#define hTestHisto_h 1

//---------------------------------------------------------------------------
//
// ClassName:   hTestHisto
//  
// Description: Singleton class to hold Emc geometry parameters.
//              User cannot access to the constructor. 
//              The pointer of the only existing object can be got via 
//              hTestHisto::GetPointer() static method. 
//              The first invokation of this static method makes 
//              the singleton object.
//
// Author:      V.Ivanchenko 27/09/00
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>
#include "G4DynamicParticle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AIDA::IHistogram1D;
class AIDA::ITuple;
class AIDA::ITree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestHisto
{

public:
  // With description

    static hTestHisto* GetPointer();

    hTestHisto();

   ~hTestHisto();
 
    void BeginOfHisto(G4int);
  // In this method histogramms are booked

    void EndOfHisto();
  // In this method bookHisto method is called in which histogramms are filled

public: // Without description

    void SetHistoName(const G4String& name) {histName = name;};
    void bookHisto();
    inline ITuple* GetNtuple() const {return ntup;};
    void SaveToTuple(const G4String&, G4double);
    void SaveToTuple(const G4String&, G4double, G4double);
    void SaveEvent();
    G4double GetTrackLength() const {return trackLength;};
    void ResetTrackLength() {trackLength = 0.0, trackAbs = true;};
    void SetTrackOutAbsorber() {trackAbs = false;};
    G4bool GetTrackInAbsorber() const {return trackAbs;};
    void AddTrackLength(G4double x)   {trackLength += x;};
    void AddEndPoint(G4double);
    void AddEnergy(G4double, G4double);
    void AddDeltaElectron(const G4DynamicParticle*);
    void AddPhoton(const G4DynamicParticle*);
    void SetVerbose(G4int val) {verbose = val;};
    G4int GetVerbose() const {return verbose;};
    void SetHistoNumber(G4int val) {nHisto = val;};
    void SetNtuple(G4bool val) {nTuple = val;};

    void SetNumberOfAbsorbers(G4int val) {NumberOfAbsorbers = val;};     
    G4int GetNumberOfAbsorbers() const {return NumberOfAbsorbers;};
    void SetAbsorberThickness(G4double val) {AbsorberThickness = val;};     
    G4double  GetAbsorberThickness() const {return AbsorberThickness;};
    void SetGap(G4double val) {gap = val;};     
    G4double  GetGap() const {return gap;};
    void SetNumAbsorbersSaved(G4int val) {nAbsSaved = val;};
    G4int GetNumAbsorbersSaved() const {return nAbsSaved;};
    void SetMaxEnergy(G4double val) {maxEnergy = val;};     
    G4double  GetMaxEnergy() const {return maxEnergy;};

private:

    hTestHisto(hTestHisto&);
    const hTestHisto& operator=(const hTestHisto&);

  // MEMBERS
    static hTestHisto* fManager;

    G4String histName;
    G4String theName;
    std::vector<AIDA::IHistogram1D*> histo;
    AIDA::ITuple* ntup;
    AIDA::ITree* tree;
    G4int nHisto;
    G4int verbose; 
    G4double zend;
    G4double zend2;
    G4double zEvt;
    G4double AbsorberThickness;
    G4double gap;
    G4int NumberOfAbsorbers;
    G4int nAbsSaved;
    G4double maxEnergy;
    G4double trackLength;
    G4bool trackAbs;
    G4bool nTuple;
};

#endif
