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
#include "g4std/fstream"
#include "g4std/iomanip"
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
  string fname1 = "dsde_1.dat";
  fout_a->open(fname1.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_b = new ofstream();
  string fname2 = "dsdtet_1.dat";
  fout_b->open(fname2.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_c = new ofstream();
  string fname3 = "dsdedtet_1.dat";
  fout_c->open(fname3.c_str(), std::ios::out|std::ios::trunc);

  ofstream* fout_a1 = new ofstream();
  string fname4 = "dsde_1.out";
  fout_a1->open(fname4.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_b1 = new ofstream();
  string fname5 = "dsdtet_1.out";
  fout_b1->open(fname5.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_c1 = new ofstream();
  string fname6 = "dsdedtet_1.out";
  fout_c1->open(fname6.c_str(), std::ios::out|std::ios::trunc);

  //there can't be lines longer than nmax characters
  const int nmax = 200;
  char line[nmax]; 
  std::string line1, line2, word1, word2, word3;
  G4bool end = true;

  G4DataVector* energy = new G4DataVector();
  int nbin = 0;
  int ibin, inum;
  int counter = 0;
  double elim = 1.07*MeV;
  double elim0= 2.0*MeV;
  double x, an, e1, e2, y1, y2, ct1, ct2, xs;
  G4DataVector* angle = new G4DataVector();
  std::vector<G4DataVector*> cs;

  // main loop 

  do {

    // read next line

    counter++;
    for(int ii = 0; ii < nmax; ii++) {line[ii] = ' ';}
    fin->getline( line, nmax);
    line1 = std::string("");
    line1 = std::string(line, nmax);
    if(1 < verbose) {
      cout << "Next line # " << counter << ": " << line1 << endl;
      cout << "First symbols= " << line[0] << line[1] << line[2] << endl;
    }  

    // analize line contence

    // end of file
    if(fin->eof()) {
      end = false;
    } else if(line[0] == 'A' && line[1] == 'D' && line[2] == 'E') {

      G4bool enddata = true;
      G4DataVector* cross = 0;
      ibin= 0;
      inum= 0;

      do {
        (*fin) >> an >> e1 >> e2 >> x >> xs;
        an *= degree;
        x *= 1000.0;      
        if(1 < verbose) {
          cout << "an= " << an/degree << " e1= " << e1 
               << " e2= " << e2 << " cross= " << x
               << " +- " << xs << " %" 
               << endl;
        }  

        if (e1 == 0.0) {
          enddata = false;
          end = false;
	}

        if((xs == 0.0 && an/degree > 10.) || !enddata) {
          // fill the rest by zero cross section
          for(int j=ibin; j<nbin; j++) {
            cross->push_back(0.0);
          }
          cs.push_back(cross);
          if(0 < verbose) {
            (*fout_c) << "#####..Result.of.parcing..####### n= " << nbin 
                      << " Angle(degree)= " << (*angle)[angle->size()-1]/degree
                      << G4endl;
            (*fout_c1) << "#####..Result.of.parcing..####### n= " << nbin 
                      << " Angle(degree)= " << (*angle)[angle->size()-1]/degree
                      << G4endl;
            for(int i=0; i<nbin; i++) {
               (*fout_c) << "e(MeV)= " << (*energy)[i]  
                         << " cross(mb/MeV/sr)= " << (*cross)[i] << endl;
               (*fout_c1) << (*cross)[i] << " ";
	    }
	    (*fout_c1) << " " << endl;
          }  
	}

	// new data
        if(xs == 0.0) {
          angle->push_back(an);
          cross = new G4DataVector();
          ibin = 0;

        } else if(end && e1 > elim) {

          ibin++;
          cross->push_back(x);
          if(an < 10.*degree) {
            energy->push_back(e1);
            nbin++;
	  }
        }      

      } while (enddata);
     
      angle->push_back(pi);
      int na = angle->size();
      G4DataVector* f1;
      G4DataVector* f2;
  
      if(0 < verbose) {
        (*fout_a) << "#####..Result.of.integration..#####.."
                  << G4endl;
        (*fout_b) << "#####..Result.of.integration..#####.. Elim(MeV)= " 
                  << elim0/MeV
                  << G4endl;
        (*fout_a1) << "#####..Result.of.integration..Energy points Nbins= " << nbin
                   << G4endl;
        (*fout_b1) << "#####..Result.of.integration..#####.. Elim(MeV)= " 
                  << elim0/MeV
                  << G4endl;

        for(int ii=0; ii<nbin; ii++) {
          (*fout_a1) << (*energy)[ii] << " ";
	}
        (*fout_a1) << G4endl;
        (*fout_a1) << "#####..Result.of.integration" << G4endl;

        for(int i=0; i<nbin; i++) {
        
          x = 0.0;
          for(int j=0; j<na-1; j++) {
            f1  = cs[j];  
            y1  = (*f1)[i];  
            ct1 = cos((*angle)[j]);
            ct2 = cos((*angle)[j+1]);
            if(j == 0) {
              f2  = cs[j+1]; 
              y2  = (*f2)[i]; 
              ct1 = 1.0;
	    } else if (j == na-2) {
              f2  = cs[j-1]; 
              y2  = (*f2)[i]; 
              ct2 = cos((*angle)[j-1]);
              y2 -= (y2 - y1)*(ct2 + 1.0)/(ct2 - ct1);
              ct2 = -1.0;
              if(y2 < 0.0) y2 = 0.0;
            } else {
              f2  = cs[j+1]; 
              y2  = (*f2)[i]; 
	    }
            x  += 0.5*(y1 + y2)*(ct1 - ct2);  
          }        
          x *= twopi;
          (*fout_a) << "e(MeV)= " << (*energy)[i] 
                    << " cross(mb/MeV)= " << x << endl;
          (*fout_a1) << x << " ";
        }

        (*fout_a1) << G4endl;

        for(int j=0; j<na-1; j++) {
          f1  = cs[j];  
          an  = cos((*angle)[j]);
          x   = 0.0;
          for(i=0; i<nbin-1; i++) {
            y1  = (*f1)[i];
            e1  = (*energy)[i];
            e2  = (*energy)[i+1];
            if(e2 < elim0) {
              y1 = 0.0;
            } else if (e1 < elim0) {
              y1 *= (e2 - elim0)/(e2 - e1);
	    }
            x += y1*(e2 - e1);
          }
          (*fout_b) << "cos(theta)= " << an 
                    << " cross(mb/sr)= " << x << endl;
          (*fout_b1) << x << " ";
        }
        (*fout_b1) << G4endl ;
      }
    }
  } while (end);
  (*fout_a) << "#####..End..#####" << G4endl;
  (*fout_b) << "#####..End..#####" << G4endl;
  (*fout_c) << "#####..End..#####" << G4endl;
}








