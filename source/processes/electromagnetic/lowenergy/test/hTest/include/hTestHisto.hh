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
#include "g4std/vector"
#include "G4DynamicParticle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class HepTupleManager;
class HepTuple;
class HepHistogram;

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
    inline HepTuple* GetNtuple() const {return ntup;};
    void SaveToTuple(const G4String&, G4double);
    void SaveToTuple(const G4String&, G4double, G4double);
    void SaveEvent();
    void AddEndPoint(G4double);
    void AddEnergy(G4double, G4double);
    void AddDeltaElectron(const G4DynamicParticle*);
    void AddPhoton(const G4DynamicParticle*);
    inline void SetVerbose(G4int val) {verbose = val;};
    inline G4int GetVerbose() const {return verbose;};
    inline void SetHistoNumber(G4int val) {nHisto = val;};

    void SetNumberOfAbsorbers(G4int val) {NumberOfAbsorbers = val;};     
    inline G4int GetNumberOfAbsorbers() const {return NumberOfAbsorbers;};
    void SetAbsorberThickness(G4double val) {AbsorberThickness = val;};     
    inline G4double  GetAbsorberThickness() const {return AbsorberThickness;};
    inline void SetNumAbsorbersSaved(G4int val) {nAbsSaved = val;};
    inline G4int GetNumAbsorbersSaved() const {return nAbsSaved;};
    void SetMaxEnergy(G4double val) {maxEnergy = val;};     
    inline G4double  GetMaxEnergy() const {return maxEnergy;};


private:

    hTestHisto(hTestHisto&);
    const hTestHisto& operator=(const hTestHisto&);

  // MEMBERS
    static hTestHisto* fManager;

    G4String histName;
    G4String theName;
    G4std::vector<HepHistogram*> histo;
    HepTupleManager* hbookManager;
    HepTuple* ntup;
    G4int nHisto;
    G4int verbose; 
    G4double zend;
    G4double zend2;
    G4double zEvt;
    G4double AbsorberThickness;
    G4int NumberOfAbsorbers;
    G4int nAbsSaved;
    G4double maxEnergy;
};

#endif
