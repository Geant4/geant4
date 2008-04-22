#include <stdio.h>

#include <iostream>
#include <fstream>

#include <string>

#include <stdlib.h>
#include <math.h>

#include "G4Incl.hh"

//#define DEBUG 1
#ifdef USEROOT // Use ROOT for data analysis

#include "TROOT.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#endif // USEROOT

using namespace std;

enum type {fullrun = 0, inclrun = 1};

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
  int type = fullrun;
  bool usingRoot = false;

  // Run type: full or cascade only
  if(strcmp("cascade", argv[1]) == 0) {
    type = inclrun;
  }

  // For ASCII output file
  ofstream out;

  G4Hazard *hazard = new G4Hazard();
  G4VarNtp *varntp = new G4VarNtp();
  G4Calincl *calincl = new G4Calincl();
  G4Ws *ws = new G4Ws();
  G4Mat *mat = new G4Mat();

  // set initial values:
  // First random seed:
  //  hazard->ial = 38035;
  hazard->ial = 979678188;

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

  G4Incl *incl = new G4Incl(hazard, calincl, ws, mat, varntp);

#ifdef USEROOT 
//   char *rootfilename = new char[2000];

//   if (rootfilename == NULL)
//   {
//     cout << "Memory allocation error." << endl;
//     return -1;
//   }

  TFile *dataFile = NULL;

  //  strcpy(rootfilename, argv[6]);
  TString rootfilename(argv[7]);
//  cout <<"filename: " << rootfilename << endl;
  if(rootfilename.Contains(".root") == 1) {
    dataFile = new TFile(rootfilename, "RECREATE");
    usingRoot = true;
  }

  // A tree into which we will put the output data
  TTree *h101 = new TTree("h101", "Data from INCL4+ABLA");

  const Int_t maxpart = 255;
  // Variables for Ntuple (TTree)
  Int_t Event;
  Int_t BulletType;
  Double_t BulletE;

  Int_t Massini, Mzini; 
  Double_t Exini;
  Double_t Pcorem, Mcorem, Pxrem, Pyrem, Pzrem;
  Int_t Mulncasc, Mulnevap,Mulntot;

  Double_t Bimpact;
  Int_t Jremn, Kfis;
  Double_t Estfis;
  Int_t Izfis, Iafis;
  Int_t Ntrack;
  Int_t baryonNumber;
  Int_t Itypcasc[maxpart], Avv[maxpart], Zvv[maxpart];
  Double_t Enerj[maxpart], Plab[maxpart], Tetlab[maxpart], Philab[maxpart];
  Double_t momX[maxpart], momY[maxpart], momZ[maxpart];

  // With these branches (variables)
  h101->Branch("Event", &Event, "Event/I");
  h101->Branch("BulletType", &BulletType, "BulletType/I");
  h101->Branch("BulletE", &BulletE, "BulletE/D");

  h101->Branch("Massini", &Massini, "Massini/I");
  h101->Branch("Mzini", &Mzini, "Mzini/I"); 
  h101->Branch("Exini", &Exini, "Exini/D");
  h101->Branch("Pcorem", &Pcorem, "Pcorem/D");
  h101->Branch("Mcorem", &Mcorem, "Mcorem/D");
  h101->Branch("Pxrem", &Pxrem, "Pxrem/D");
  h101->Branch("Pyrem", &Pyrem, "Pyrem/D");
  h101->Branch("Pzrem", &Pzrem, "Pzrem/D");
  h101->Branch("Mulncasc", &Mulncasc, "Mulncasc/I");
  h101->Branch("Mulnevap", &Mulnevap, "Mulnevap/I");
  h101->Branch("Mulntot", &Mulntot, "Mulntot/I");

  h101->Branch("Bimpact", &Bimpact, "Bimpact/D");
  h101->Branch("Jremn", &Jremn, "Jremn/I");
  h101->Branch("Kfis", &Kfis, "Kfis/I");
  h101->Branch("Estfis", &Estfis, "Estfis/D");
  h101->Branch("Izfis", &Izfis, "Izfis/I");
  h101->Branch("Iafis", &Iafis, "Iafis/I");

  h101->Branch("baryonNumber", &baryonNumber, "baryonNumber/I");
  h101->Branch("Ntrack", &Ntrack, "Ntrack/I");
  h101->Branch("Itypcasc", Itypcasc, "Itypcasc[Ntrack]/I");
  h101->Branch("Avv", Avv, "Avv[Ntrack]/I"); 
  h101->Branch("Zvv", Zvv, "Zvv[Ntrack]/I");
  h101->Branch("Enerj", Enerj, "Enerj[Ntrack]/D");
  h101->Branch("Plab", Plab, "Plab[Ntrack]/D");
  h101->Branch("Tetlab", Tetlab, "Tetlab[Ntrack]/D");
  h101->Branch("Philab", Philab, "Philab[Ntrack]/D");
  h101->Branch("momX", momX, "momX[Ntrack]/D");
  h101->Branch("momY", momY, "momY[Ntrack]/D");
  h101->Branch("momZ", momZ, "momZ[Ntrack]/D");
  
#endif // USEROOT

  if (filename == NULL)
  {
    cout << "Memory allocation error." << endl;
    return -1;
  }

  if(argc < 8) {
    cout <<"Usage: inclStandAlone targetA targetZ bulletType bulletEnergy events outputfile" << endl;
    return -1;
  }

  if(!usingRoot) {
    strcpy (filename, argv[7]);
    // Open the output file:
    out.open(filename);
  }

  // Initialize FINPUT:
  for(int i = 0; i < 15; i++) {
    calincl->f[i] = 0.0;
  }
  // End of Initialization

  // Input parameters:
  // targetA, targetZ, bulletType, bulletE, filename
  // Usage: cppinterfacetest targetA targetZ bulletType bulletE filename

  // Set input parameters:

  // Target nucleus:
  // Mass number:
  //FINPUT(1)
  calincl->f[0] = atof(argv[2]);
  // Charge number:
  // FINPUT(2)
  calincl->f[1] = atof(argv[3]);

  // Bullet:
  // Bullet type (1.00 = proton):
  // FINPUT(7)
  calincl->f[6] = atof(argv[4]);
  // Bullet energy:
  // FINPUT(3)
  calincl->f[2] = atof(argv[5]);

  // Run settings:
  // Time scaling:
  // FINPUT(6)
  calincl->f[5] = 1.0;

  // Evaporation model 1 = KHS, 2 = GEM:
  //  localvars_.choice_evap = 1;

  // Other random seeds:
  //   int seeds[19] = {21033,17563,27563,33657,43657,56375,66375,77365,87365, 
  // 		   12345,22345,32345,42345,52345,62345,72345,82345,34567,
  // 		   47059};
  //   for(int seedIndex = 0; seedIndex < 19; seedIndex++) {
  //     hazard->igraine[seedIndex] = seeds[seedIndex];
  //   }

  //   cout <<"Seeds at the beginning of the simulation:" << endl;
  //   cout <<"ial: " << hazard_.ial << endl;
  //   for(int seedIndex = 0; seedIndex < 19; seedIndex++) {
  //     cout <<"IY(" << seedIndex <<") = " << hazard_.IY[seedIndex] << endl;
  //   }
  //   cout << endl;


  // Nuclear potential:
  // FINPUT(5)
  calincl->f[4] = 45.0;

  // NOSURF:
  ws->nosurf = -2;

  // XFOISA
  ws->xfoisa = 8;

  // NPAULSTR
  ws->npaulstr = 0;

  // Message output file:
  //  sprintf(localvars_.stringout, "inclabla.out");

  // Data path:
  //  sprintf(localvars_.RACINE, "./dapnia/dapx4");

  // Events (only for compatibility)
  // Deprecated
  calincl->icoup = 1;

  // Events:
  //  int totalevents = 100000;
  int totalevents = atoi(argv[6]);
  //  debugval_.allevents = totalevents;

  int particleI = 0;
  
  // End of input parameters

  cout << "Outpufile: " << argv[7] << endl;
  cout << "Bullet: " << endl;
  cout << "Type: " << calincl->f[6] << endl;
  cout << "energy: " << calincl->f[2] << " mev " << endl;
  cout << "target: " << endl;
  cout << "a: " << calincl->f[0] << endl;
  cout << "z: " << calincl->f[1] << endl;
  cout << "events: " << argv[6] << endl;
  cout << endl;

  cout << "Running..." << endl;

  mat->nbmat = 1;
  mat->amat[0] = int(calincl->f[0]);
  mat->zmat[0] = int(calincl->f[1]);

  cout <<"Initializing INCL" << endl;
  incl->initIncl(false);
  cout <<"Initialization complete." << endl;

  for(int n = 1; n <= totalevents; n++) {
    if(n == 1) {
      doinit = 1;
    } 
    else {
      doinit = 0;
    }

    if(type == fullrun) {
      incl->processEventInclAbla(n);
    } else if(type == inclrun) {
      incl->processEventIncl();
    } else {
      cout <<"Unknown run type." << endl;
      exit(1);
    }
    // The coulomb transparency:
    if(n == 1) {
      //      totalevents = totalevents - debugval_.ntranscoul;
      totalevents = totalevents - 1;
    }

    //   *      1   * I*4  *             * VARNTP  * massini	     A of the remnant
    // *      2   * I*4  *             * VARNTP  * MZINI	     Z    "        "
    // *      3   * R*4  *             * VARNTP  * EXINI	     Excit energy " "
    // *      4   * I*4  *             * VARNTP  * MULNCASC	     Cascade n multip.
    // *      5   * I*4  *             * VARNTP  * MULNEVAP	     Evapo   "      "
    // *      6   * I*4  *             * VARNTP  * MULNTOT	     Total   "      "
    // *      7   * R*4  *             * VARNTP  * BIMPACT	     Impact parameter
    // *      8   * I*4  *             * VARNTP  * JREMN	     Remnant Intrinsic Spin
    // *      9   * I*4  *             * VARNTP  * KFIS	     Fission 1/0=Y/N
    // *     10   * R*4  *             * VARNTP  * ESTFIS		Excit energy at fis
    // *     11   * I*4  *             * VARNTP  * IZFIS		Z of fiss nucleus
    // *     12   * I*4  *             * VARNTP  * IAFIS		A of "          "
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
      //      std::cout <<"Filling event number: " << n << std::endl;
      double momXsum = 0.0, momYsum = 0.0, momZsum = 0.0;
      for(particleI = 0; particleI < varntp->ntrack; particleI++) {
#ifdef DEBUG
	  momX[particleI] = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Cos(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momY[particleI] = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Sin(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momZ[particleI] = varntp->plab[particleI]*TMath::Cos(varntp->tetlab[particleI]*TMath::Pi()/180.0);
	  momXsum += momX[particleI];
	  momYsum += momY[particleI];
	  momZsum += momZ[particleI];
	  G4cout <<"A = " << varntp->avv[particleI] << " Z = " << varntp->zvv[particleI] << " mom = (" << momX[particleI] << ", " << momY[particleI] << ", " << momZ[particleI] <<")" << G4endl;
#endif
	if(!usingRoot) {
	  out << n << " " << calincl->f[6]  << " " << calincl->f[2] << " ";
	  out << varntp->massini << " " << varntp->mzini << " ";
	  out << varntp->exini << " " << varntp->mulncasc << " " << varntp->mulnevap << " " << varntp->mulntot << " ";
	  out << varntp->bimpact << " " << varntp->jremn << " " << varntp->kfis << " " << varntp->estfis << " ";
	  out << varntp->izfis << " " << varntp->iafis << " " << varntp->ntrack << " ";
	  out << varntp->itypcasc[particleI] << " ";
	  out << varntp->avv[particleI] << " " << varntp->zvv[particleI] << " " << varntp->enerj[particleI] << " ";
	  out << varntp->plab[particleI] << " " << varntp->tetlab[particleI] << " " << varntp->philab[particleI] << endl;
	}
#ifdef USEROOT
	if(usingRoot) {
	  //	  std::cout <<"Filling particle number: " << particleI << std::endl;
	  Event = n;
	  BulletType = (Int_t)calincl->f[6]; 
	  BulletE = calincl->f[2];
	  Massini = int(varntp->massini);
	  Mzini = int(varntp->mzini);
	  Exini = varntp->exini;
	  Pcorem = varntp->pcorem;
	  Mcorem = varntp->mcorem;
	  Pxrem = varntp->pxrem;
	  Pyrem = varntp->pyrem;
	  Pzrem = varntp->pzrem;
	  Mulncasc = varntp->mulncasc;
	  Mulnevap = varntp->mulnevap;
	  Mulntot = varntp->mulntot;
	  Bimpact = varntp->bimpact; 
	  Jremn = varntp->jremn; 
	  Kfis = varntp->kfis;
	  Estfis = varntp->estfis;
	  Izfis = varntp->izfis; 
	  Iafis = varntp->iafis; 
	  baryonNumber = varntp->getTotalBaryonNumber();
	  Ntrack = varntp->ntrack;
	  Itypcasc[particleI] = varntp->itypcasc[particleI];
	  Avv[particleI] = varntp->avv[particleI];
	  Zvv[particleI] = varntp->zvv[particleI];
	  Enerj[particleI] = varntp->enerj[particleI];
	  Plab[particleI] = varntp->plab[particleI];
	  Tetlab[particleI] = varntp->tetlab[particleI];
	  Philab[particleI] = varntp->philab[particleI];
	  momX[particleI] = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Cos(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momY[particleI] = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Sin(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momZ[particleI] = varntp->plab[particleI]*TMath::Cos(varntp->tetlab[particleI]*TMath::Pi()/180.0);
	  
	  //	  std::cout <<"Particle " << particleI << " of " << Ntrack - 1<<" filled." << std::endl;
	}
#endif // USEROOT
      }
#ifdef DEBUG
      G4cout <<"-------------------------------------------------" << G4endl;
      G4cout <<" mom = (" << momXsum << ", " << momYsum << ", " << momZsum << ")" << G4endl;
      G4cout <<"Total momentum: " << varntp->getMomentumSum() << G4endl;
#endif
#ifdef USEROOT 
      if(usingRoot) {
	//	std::cout <<"Filling..." << std::endl;
	h101->Fill();
	//	std::cout <<"Filling complete." << std::endl;
      }
#endif // USEROOT
      varntp->ntrack = 0;
    }
    else {
      if(varntp->ntrack == -2) {
	n = n - 1;
      }
      else {
	transparent++;
      }
    }
  }
  
  if(!usingRoot) {
    out.close();
  }

  sprintf(summaryFilename, "%s.runSummary", (char*)argv[7]);

  ofstream summaryFile;
  summaryFile.open(summaryFilename);
  
  summaryFile << "INCL4+ABLA C++-Fortran hybrid thin-target calculation:" << endl;
  summaryFile << endl;
  summaryFile << "Run setup:" << endl;
  summaryFile << "Bullet: " << endl;
  summaryFile << "\t Type: " << calincl->f[6] << endl;
  summaryFile << "\t Energy: " << calincl->f[2] << " MeV " << endl;
  summaryFile << "Target: " << endl;
  summaryFile << "\t A: " << calincl->f[0] << endl;
  summaryFile << "\t Z: " << calincl->f[1] << endl;
  summaryFile << "Events: " << argv[6] << endl;
  summaryFile << endl;
  if(!usingRoot) {
    summaryFile << "Calculation output in ASCII file: " << filename << endl;
  }
#ifdef USEROOT
  else {
    summaryFile << "Calculation output in ROOT file: " << rootfilename << endl;
  }
#endif // USEROOT
  summaryFile << endl;

//         f_cross_sect=(icoup-ntrans)
// 	f_cross_sect=f_cross_sect/(icoup+ntrans_coul)
  fCrossSection = (((double)totalevents) - ((double)transparent))/((double) totalevents + 0); //(double) debugval_.ntranscoul);


//  geomCrossSection = fCrossSection * 31.4159 * ws_.BMAX * ws_.BMAX;
  geomCrossSection = 31.4159*(ws->bmax)*(ws->bmax);

  summaryFile << "Output information: " << endl;
  summaryFile << "Total number of events: " << endl;
  summaryFile << "\t Asked by the user: " << (totalevents + 0 ) << endl; //debugval_.ntranscoul) << endl;
  summaryFile << "\t Transparent: " << transparent << endl;
  summaryFile << "\t Proper cascades: " << (totalevents - transparent) << endl;
  summaryFile << "Maximum impact parameter: " << ws->bmax << endl;
  summaryFile << "Geometrical cross-section: " << geomCrossSection << " mb" << endl;
  summaryFile << "Total reaction cross section: " << fCrossSection*31.4159*pow(ws->bmax,2) << " mb" << endl;
  summaryFile << "---------------------" << G4endl;
  summaryFile << "Normalization factor = " << geomCrossSection/(double(totalevents) - double(transparent)) << endl;
  summaryFile.close();

#ifdef USEROOT
  if(dataFile != NULL) {
    std::cout <<"Writing data..." << std::endl;
    h101->Write();
    std::cout <<"Closing ROOT file..." << std::endl;
    dataFile->Close();
  }
#endif // USEROOT
  
  cout << "Simulation run done." << endl;
  return 0;
}

