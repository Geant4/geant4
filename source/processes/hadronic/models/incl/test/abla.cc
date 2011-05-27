#include <stdio.h>

#include <iostream>
#include <fstream>

#include <string>

#include <stdlib.h>
#include <math.h>

#include "G4Incl.hh"

#include "G4NucleiProperties.hh"

#ifdef USEROOT // Use ROOT for data analysis

#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#endif // USEROOT

using namespace std;

int main(int argc, char *argv[])
{
  //  const char outputfile[100] = argv[5];

  char *filename = new char [2000];
  char *summaryFilename = new char [2000];

  // Summary and diagnostics data:
  double fCrossSection;
  double geomCrossSection;

  int transparent = 0;

  int doinit = 1;

  bool usingRoot = true;

  // For ASCII output file
  ofstream out;

  G4Hazard *hazard = new G4Hazard();
  G4VarNtp *varntp = new G4VarNtp();
  G4Volant *volant = new G4Volant();
  
  // set initial values:
  // First random seed:
  // (Premiere graine)
  //  hazard->ial = 38035;
  hazard->ial = 979678188;
  //  hazard->ial = 387553;
  // other seeds:
  hazard->igraine[0] = 3997;
  hazard->igraine[1] = 15573;
  hazard->igraine[2] = 9971;
  hazard->igraine[3] = 9821; 
  hazard->igraine[4] = 99233; 
  hazard->igraine[5] = 11167; 
  hazard->igraine[6] = 12399;
  hazard->igraine[7] = 11321; 
  hazard->igraine[8] = 9825;
  hazard->igraine[9] = 2587; 
  hazard->igraine[10] = 1775;
  hazard->igraine[11] = 56799; 
  hazard->igraine[12] = 1156;
  //  hazard->igraine[13] = 11207;
  hazard->igraine[13] = 38957; 
  hazard->igraine[14] = 35779; 
  hazard->igraine[15] = 10055; 
  hazard->igraine[16] = 76533; 
  hazard->igraine[17] = 33759;
  hazard->igraine[18] = 13227;

  //  G4Incl *incl = new G4Incl(hazard, calincl, ws, mat, varntp);
  G4Abla *abla = new G4Abla(hazard, volant, varntp);
  
  // Open the output file:
  out.open("events.out");

  // Input:
  G4double nucleusA = atof(argv[1]);
  G4double nucleusZ = atof(argv[2]);
  G4double excitationE = atof(argv[3]);
  double angularMom = atof(argv[4]);
  G4int numberOfEvents = atoi(argv[5]);
  G4double mass = G4NucleiProperties::GetNuclearMass(nucleusA,nucleusZ)/MeV;

  G4double recoilE = 0;
  G4double nucMomX = 0;
  G4double nucMomY = 0;
  G4double nucMomZ = 0;

  cout << "Outpufile: events.out " << endl;
  cout << "Nucleus A: " << nucleusA << endl;
  cout << "Nucleus Z: " << nucleusZ << endl;
  cout << "Atomic mass (from Geant4): " << mass << endl;  
  cout << "Excitation energy: " << excitationE << endl;
  cout << endl;

  Double_t Avv, Zvv, Enerj, Plab, Tetlab, Philab;

  TFile *f = new TFile("abla.root", "RECREATE");
  //  TNtuple *data = new TNtuple("data", "ABLA output", "Avv/F:Zvv/F:Enerj/F:Plab/F:Tetlab/F:Philab/F");
  TNtuple *data = new TNtuple("data", "ABLA output", "Avv:Zvv:Enerj:Plab:Tetlab:Philab");
  abla->initEvapora();

  double momX, momY, momZ;
  double momXsum, momYsum, momZsum;

  for(G4int event = 1; event <= numberOfEvents; event++) {
    momX = 0.0;
    momY = 0.0;
    momZ = 0.0;
    momXsum = 0.0;
    momYsum = 0.0;
    momZsum = 0.0;
    abla->breakItUp(nucleusA, nucleusZ, mass, excitationE, angularMom,
		    recoilE, nucMomX, nucMomY, nucMomZ, event);
    // *     13   * I*4  *[0,250]      * VARNTP  * NTRACK		Number of particles
    // *     14   * I*4  *             * VARNTP  * ITYP(NTRACK)	emitted in cascade (0)
    // *                                                                       or evapo   (1)
    // *     15   * I*4  *             * VARNTP  * AVV(NTRACK)		A (-1 for pions)
    // *     16   * I*4  *             * VARNTP  * ZVV(NTRACK)		Z
    // *     17   * R*4  *             * VARNTP  * ENERJ(NTRACK)		kinetic energy
    // *     18   * R*4  *             * VARNTP  * PLAB(NTRACK)		momentum
    // *     19   * R*4  *             * VARNTP  * TETLAB(NTRACK)		Theta (deg)
    // *     20   * R*4  *             * VARNTP  * PHILAB(NTRACK)		Phi   (deg)
    if(varntp->ntrack > 0) {
      for(G4int particleI = 0; particleI < varntp->ntrack; particleI++) {
 	momX = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Cos(varntp->philab[particleI]*TMath::Pi()/180.0);
 	momY = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Sin(varntp->philab[particleI]*TMath::Pi()/180.0);
 	momZ = varntp->plab[particleI]*TMath::Cos(varntp->tetlab[particleI]*TMath::Pi()/180.0);
// 	momX = varntp->plab[particleI]*TMath::Sin(varntp->philab[particleI]*TMath::Pi()/180.0)*TMath::Cos(varntp->tetlab[particleI]*TMath::Pi()/180.0);
// 	momY = varntp->plab[particleI]*TMath::Sin(varntp->philab[particleI]*TMath::Pi()/180.0)*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0);
// 	momZ = varntp->plab[particleI]*TMath::Cos(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momXsum += momX;
	  momYsum += momY;
	  momZsum += momZ;

	  G4cout <<"A = " << varntp->avv[particleI] << " Z = " << varntp->zvv[particleI] << " mom = (" << momX << ", " << momY << ", " << momZ <<")" << G4endl;
	out << event << " " << varntp->avv[particleI] << " " 
	    << varntp->zvv[particleI] << " " << varntp->enerj[particleI] << " "
	    << varntp->plab[particleI] << " " << varntp->tetlab[particleI] << " " 
	    << varntp->philab[particleI] << endl;
	Avv = varntp->zvv[particleI];
	Zvv = varntp->avv[particleI];
	Enerj = varntp->enerj[particleI];
	Plab= varntp->plab[particleI];
	Tetlab= varntp->tetlab[particleI];
	Philab = varntp->philab[particleI];
	data->Fill(Avv, Zvv, Enerj, Plab, Tetlab, Philab);
      }
      G4cout <<"-------------------------------------------------" << G4endl;
      G4cout <<" mom = (" << momXsum << ", " << momYsum << ", " << momZsum << ")" << G4endl;
      G4cout <<"Total momentum: " << varntp->getMomentumSum() << G4endl;
      varntp->ntrack = 0;
    }
    else {
      transparent++;
    }
    varntp->ntrack = -1;
  }
  
  out.close();
  data->Write();
  f->Close();

  cout << "Done." << endl;
  return 0;
}

