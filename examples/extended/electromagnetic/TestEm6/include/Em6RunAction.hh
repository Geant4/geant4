// Em6RunAction.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6RunAction_h
#define Em6RunAction_h 1

#include "G4UserRunAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "g4std/vector"

typedef  G4std::vector<G4double> MyVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6DetectorConstruction;
class Em6PrimaryGeneratorAction;

class G4Run;

#ifndef G4NOHIST
 class HepTupleManager;
 class HepHistogram;
#endif

class Em6RunAction : public G4UserRunAction
{
  public:

    Em6RunAction(Em6DetectorConstruction*, Em6PrimaryGeneratorAction*);
   ~Em6RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);

    inline void initializePerEvent();
           void fillPerEvent();
    inline void fillPerTrack(G4double,G4double);
    inline void fillPerStep (G4double,G4int);
    inline void particleFlux(G4ParticleDefinition*,G4int);
    inline void fillGamma(G4int GammaID,G4double Egam,G4ThreeVector Pgam)
	  { if(this->Egam==0) // store only for the first gamma in the event with -> mu+mu-
        { this->GammaID=GammaID;
          this->Egam=Egam;
	      this->Pgam=Pgam;
	    }
      }
    inline void fillMuPlus(G4double EMuPlus,G4ThreeVector PMuPlus)
	  { this->EMuPlus=EMuPlus;
        this->PMuPlus=PMuPlus;
	    xPlus=EMuPlus/Egam;
	    thetaPlus=Pgam.angle(PMuPlus); // angle between gamma and mu+
      }
    inline void fillMuMinus(G4double EMuMinus,G4ThreeVector PMuMinus)
	  { this->EMuMinus=EMuMinus;
        this->PMuMinus=PMuMinus;
        thetaMinus=Pgam.angle(PMuMinus); // angle between gamma and mu-
      }
    inline void IncrementGammaMuMuCounter() { NGamMuMu++; };
    inline G4int GetGammaID() { return GammaID; }

  private:

    void bookHisto();
    void cleanHisto();

  private:

    Em6DetectorConstruction*   Em6Det;
    Em6PrimaryGeneratorAction* Em6Kin;

    G4int nLbin;
    MyVector dEdL;
    MyVector sumELongit;
    MyVector sumE2Longit;
    MyVector sumELongitCumul;
    MyVector sumE2LongitCumul;

    MyVector gammaFlux;
    MyVector electronFlux;
    MyVector positronFlux;

    G4double ChargTrLength;
    G4double sumChargTrLength;
    G4double sum2ChargTrLength;

    G4double NeutrTrLength;
    G4double sumNeutrTrLength;
    G4double sum2NeutrTrLength;

    // info on the gamma->mu+mu- process  (if several the first in the event)
    G4ThreeVector Pgam,PMuPlus,PMuMinus; // photon and muon momentum vectors
    G4int GammaID,NGamMuMu;
    G4double Egam,EMuPlus,EMuMinus, //photon and muon energies
	         xPlus,                 // Emu+/Egam energye
	         thetaPlus,thetaMinus;  //photon energy in last gamma->mu+mu- process

#ifndef G4NOHIST
    HepTupleManager* hbookManager;
    HepHistogram *histo1, *histo2, *histo3;
    HepHistogram *histo4, *histo5, *histo6;
    HepHistogram *histo7, *histo8, *histo9;
    HepHistogram *hist10, *hist11, *hist12;
    HepHistogram *hist13, *hist14;
#endif

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Em6RunAction::initializePerEvent()
{
  //initialize arrays of energy deposit per bin
  for (G4int i=0; i<nLbin; i++)
     { dEdL[i] = 0.; }

  //initialize tracklength
    ChargTrLength = NeutrTrLength = 0.;

    Egam=EMuPlus=EMuMinus=xPlus=thetaPlus=thetaMinus=0; // init gamma->mu+mu- info
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Em6RunAction::fillPerTrack(G4double charge, G4double trkLength)
{
  if (charge != 0.) ChargTrLength += trkLength;
  else              NeutrTrLength += trkLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Em6RunAction::fillPerStep(G4double dEstep, G4int Lbin)
{
  dEdL[Lbin] += dEstep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void Em6RunAction::particleFlux(G4ParticleDefinition* particle, G4int Lplan)
{
       if (particle == G4Gamma::Gamma())          gammaFlux[Lplan]++;
  else if (particle == G4Electron::Electron()) electronFlux[Lplan]++;
  else if (particle == G4Positron::Positron()) positronFlux[Lplan]++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif

