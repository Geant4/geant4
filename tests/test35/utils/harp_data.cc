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
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     Test30
//
//      Author:        V.Ivanchenko 
// 
//      Creation date: 19 February 2006
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4ios.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

void Print(std::ofstream* fout, std::string s, std::string s1, std::string s2,
	   int n, double* x) {
  (*fout) << "**  " << s << std::endl;
  (*fout) << "ve/create " << s1 << s2 << "(" << n << ") R ";
  for(int i=0; i<n; i++) {(*fout) << " " << x[i];}
  (*fout) << std::endl;
}

void Print1(std::ofstream* fout, std::string s, std::string s1, std::string s2,
	   int n, double* x) {
  (*fout) << "**  " << s << std::endl;
  (*fout) << "ve/create " << s1 << s2 << "(" << 2*n+1 << ") R ";
  if(s2 == "x") (*fout) << " 0  ";
  for(int i=0; i<n; i++) {(*fout) << " " << x[i] << " " << x[i];}
  if(s2 != "x") (*fout) << " 0  ";
  (*fout) << std::endl;
}

int main(int argc, char** argv)
{
  int verbose = 1;

  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  std::ifstream* fin = new std::ifstream();
  std::string fname = argv[1];
  std::string finName = fname + ".txt";
  fin->open(finName.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

  std::ofstream* fout = new std::ofstream();
  std::string foutName = fname + ".dat";
  fout->open(foutName.c_str(), std::ios::out|std::ios::trunc);

  //there can't be lines longer than nmax characters
  const int nmax = 200;
  char line[nmax]; 
  std::string line1, line2, word1, word2, word3;
  G4bool end = false;

  double mom[50];
  double ang[50];
  double cosa[50];
  double m[50];
  double a[50];
  double sm[50];
  double sa[50];
  double cross[50][50];
  double rms[50][50];
  double cross_m[50];
  double rms_m[50];
  double cross_a[50];
  double rms_a[50];
  int nmom = 0;
  int nang = 0;
  int counter = 0;

  // main loop 

  do {

    // read next line
    counter++;
    for(int ii = 0; ii < nmax; ii++) {line[ii] = ' ';}
    fin->getline( line, nmax);
    line1 = std::string("");
    line1 = std::string(line, nmax);
    if(1 < verbose) {
      std::cout << "Next line # " << counter << ": " << line1 << std::endl;
      std::cout << "First symbols= " 
		<< line[0] << line[1] << line[2] << std::endl;
    }  

    // analize line contence

    // end of file
    if(line[0] == '#' && line[1] == 'e' && line[2] == 'n') {
      end = true;
    } else if(line[0] == '#' && line[1] == 'a' && line[2] == 'n') {

      (*fin) >> nang;
      if(1 < verbose) std::cout << "### Vector of Angles (radian)" << std::endl;
      for(int i=0; i<nang; i++) {
        (*fin) >> ang[i];
        cosa[i] = std::cos(ang[i]);
	if(1 < verbose) std::cout << "  " << ang[i];
      }
      if(1 < verbose) std::cout << std::endl;
      nang--;
    } else if(line[0] == '#' && line[1] == 'm' && line[2] == 'o') {

      (*fin) >> nmom;
      if(1 < verbose) std::cout << "### Vector of Momentum (Gev/c)" << std::endl;
      for(int i=0; i<nmom; i++) {
        (*fin) >> mom[i];
	if(1 < verbose) std::cout << "  " << mom[i];
      }
      if(1 < verbose) std::cout << std::endl;
      nmom--;
    } else if(line[0] == '#' && line[1] == 'c' && line[2] == 'r') {

      if(1 < verbose) std::cout << "### Table of Cross Section (mb/GeV/c sr))" << std::endl;
      for(int i=0; i<nang; i++) {
	for(int j=0; j<nmom; j++) {
	  (*fin) >> cross[i][j] >> rms[i][j];
	  if(1 < verbose) std::cout << "  " << cross[i][j] << " +- " << rms[i][j] << std::endl;
	}
      }
    }
  } while (!end);

  double r, rs, cr, scr;
  int i, j;
  for(i=0; i<nang; i++) {
    a[i] = 0.5*(ang[i+1] + ang[i]);
    sa[i]= a[i] - ang[i];
  }
  for(i=0; i<nmom; i++) {
    m[i] = 0.5*(mom[i+1] + mom[i]);
    sm[i]= m[i] - mom[i];
    rs = 0.0;
    cr = 0.0;
    for(j=0; j<nang; j++) {
      double f = (cosa[j] - cosa[j+1])*twopi;
      cr += cross[j][i]*f;
      r   =  rms[j][i]*f;
      rs += r*r;
    }
    cross_m[i] = cr;
    scr = cr*0.058;
    rms_m[i] = std::sqrt(rs + scr*scr);
  }
  for(i=0; i<nang; i++) {
    rs = 0.0;
    cr = 0.0;
    for(j=0; j<nmom; j++) {
      cr += cross[i][j]*2.0*sm[j];
      r   =  rms[i][j]*2.0*sm[j];
      rs += r*r;
    }
    cross_a[i] = cr;
    scr = cr*0.058;
    rms_a[i] = std::sqrt(rs + scr*scr);
  }
  Print(fout, "Angles (radians)","an","",nang,a);
  Print(fout, "RMS Angles (radians)","rmsa","",nang,sa);
  Print(fout, "Cross Section dSig/do (mb/sr)","cr_a","",nang,cross_a);
  Print(fout, "RMS Cross Section dSig/do (mb/sr)","rmscr_a","",nang,rms_a);

  Print(fout, "Momentum (Gev/c)","p","",nmom,m);
  Print(fout, "RMS Momentum (Gev/c)","rmsp","",nmom,sm);
  Print(fout, "Cross Section dSig/dp (mb/GeV/c)","cr_p","",nmom,cross_m);
  Print(fout, "RMS Cross Section dSig/dp (mb/GeV/c)","rmscr_p","",nmom,rms_m);

  std::string s[10] = {"0","1","2","3","4","5","6","7","8","9"};
  for(i=0; i<nang; i++) {
    (*fout) << "#####.. Theta= " 
	    << ang[i] << " - " << ang[i+1] << std::endl;
    Print(fout, "Cross Section dSig/dpdo (mb/(Gev/c sr))","c",s[i],nmom, &(cross[i][0]));
    Print(fout, "RMS Cross Section dSig/dpdo (mb/(Gev/c sr))","rmsc",s[i],nmom,&(rms[i][0]));
  }

  // ds/dp(mb/GeV) for pi+ with theta cut 0.03-0.21 radian LEPAR
  double s1[9]={162.877,  82.677,  62.0968,  49.4992,  38.6704,  30.4079,  22.4535,  13.5883,  5.68604};
  // ds/dp(mb/GeV) for pi+ with theta cut 0.03-0.21 radian Binary
  double s2[9]={14.1956,  26.795,  30.814,  31.8324,  30.3898,  28.1858,  24.1517,  17.0263,  8.40627};
  // ds/dp(mb/GeV) for pi+ with theta cut 0.03-0.21 radian Bertini
  double s3[9]={19.8601,  23.6435,  24.9681,  23.5527,  18.2391,  14.3827,  11.4072,  9.06188,  7.28947};
  // ds/dp(mb/GeV) for pi+ with theta cut 0.03-0.21 radian QGSP
  double s4[9]={33.4052,  66.2793,  71.5502,  60.3934,  46.7466,  34.3234,  22.852,  11.7641,  4.5756};
  // ds/dp(mb/GeV) for pi+ with theta cut 0.03-0.21 radian QGSC
  double s5[9]={41.0835,  72.1916,  74.4173,  61.7734,  47.2698,  34.7077,  22.9362,  11.949,  4.64495};
  Print1(fout, "p  ","p","x",9,mom);
  Print1(fout, "LEPAR ","s1","y",9,s1);
  Print1(fout, "Binary ","s2","y",9,s2);
  Print1(fout, "Bertini ","s3","y",9,s3);
  Print1(fout, "QGSP ","s4","y",9,s4);
  Print1(fout, "QGSC ","s5","y",9,s5);

  (*fout) << "**..End..**" << std::endl;
}








