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
//
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
//      Creation date: 12 March 2002
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <string>
#include <iostream.h>
#include <stdlib.h>
#include <strstream.h>


int main(int argc, char** argv)
{
  int verbose = 1;

  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  ifstream* fin = new ifstream();
  string fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }
  
  ofstream* fout_a = new ofstream();
  string fname1 = "dsde.out";
  fout_a->open(fname1.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_b = new ofstream();
  string fname2 = "dsdedtet.out";
  fout_b->open(fname2.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_c = new ofstream();
  string fname3 = "dsdedtet.dat";
  fout_c->open(fname3.c_str(), std::ios::out|std::ios::trunc);

  (*fout_c) << "#####..Result.of.parcing..#######.. "
            << " Angle(degree)= 0" << G4endl;


  //there can't be lines longer than nmax characters
  const int nmax = 100;
  char line[nmax];
  std::string line1, line2, word1, word2, word3;
  G4bool end = true;
  double ebeam = 0.0;
  std::vector<double> energy;
  std::vector<double> cross;
  std::vector<double> err;

  // main loop
  int counter = 0;

  do {
    counter++;
    for(int ii = 0; ii < nmax; ii++) {line[ii] = ' ';}
    fin->getline( line, nmax);
    line1 = std::string("");
    line1 = std::string(line, nmax);
    if(2 < verbose) {
      cout << "Next line # " << counter << ": " << line1 << endl;
      cout << "First symbols= " << line[0] << line[1] << line[2] << endl;
    }

    // analize line contence

    // end of file
    if(fin->eof() || line[0] == 'E' && line[1] == 'N' && line[2] == 'D') {
      end = false;
    } else if(line[0] == 'M' && line[1] == 'E' && line[2] == 'V') {
      (*fin) >> ebeam;
    } else if(line[0] == 'A' && line[1] == 'D' && line[2] == 'E') {

      double an, e0, e, s, s1;
      int i, nbin;
      (*fin) >> an;
      (*fin) >> word1 >> word2 >> nbin;
      (*fin) >> word1 >> word2 >> word3 >> line1;
      (*fin) >> word1 >> word2 >> word3 >> line1;
      energy.clear();
      cross.clear();
      err.clear();
      if(an > 0.0) {
        (*fout_b) << "#####..Ebeam(MeV) = " << ebeam/MeV << " Theta(deg)= " << an << G4endl;
        (*fout_b) << "ve/create X(" << nbin << ") R ";
        (*fout_c) << "#####..Ebeam(MeV) = " << ebeam/MeV << " Theta(deg)= " << an << G4endl;
      } else {
        (*fout_a) << "#####..Ebeam(MeV) = " << ebeam/MeV << " ds/dE " << G4endl;
        (*fout_a) << "ve/create X(" << nbin << ") R ";
        (*fout_c) << "#####..Ebeam(MeV) = " << ebeam/MeV << " ds/dE " << G4endl;
      }
      for(i=0; i<nbin; i++) {
        (*fin) >> e >> s >> e0 >> s1;
        if(1 < verbose) {
          cout << "e= " << e << " cross= " << s << " +- " << s1 << endl;
        }
        energy.push_back(e);
	cross.push_back(s);
	err.push_back(s1);
        (*fout_c) << i << ".   e(MeV)= " << e;
	if(an> 0.0) {
	  (*fout_c) << "  cross(mb/sr/MeV)= ";
          (*fout_b) << e << " ";
	} else {
	  (*fout_c) << "  cross(mb/MeV)= ";
          (*fout_a) << e << " ";
	}
	(*fout_c) << s << " +- " << s1 << endl;
      }
      if(an > 0.0) {
        (*fout_b) << endl;
        (*fout_b) << "ve/create Y(" << nbin << ") R ";
      } else {
        (*fout_a) << endl;
        (*fout_a) << "ve/create Y(" << nbin << ") R ";
      }
      for(i=0; i<nbin; i++) {
	if(an> 0.0) {
          (*fout_b) << cross[i] << " ";
	} else {
          (*fout_a) << cross[i] << " ";
	}
      }
       if(an > 0.0) {
        (*fout_b) << endl;
        (*fout_b) << "ve/create Z(" << nbin << ") R ";
      } else {
        (*fout_a) << endl;
        (*fout_a) << "ve/create Z(" << nbin << ") R ";
      }
      for(i=0; i<nbin; i++) {
	if(an> 0.0) {
          (*fout_b) << err[i] << " ";
	} else {
          (*fout_a) << err[i] << " ";
	}
      }
      (*fout_a) << endl;
      (*fout_b) << endl;
    }
  } while (end);
  (*fout_a) << "#####..End..#####" << G4endl;
  (*fout_b) << "#####..End..#####" << G4endl;
  (*fout_c) << "#####..End..#####" << G4endl;
}






