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
//      File name:     eedl
//                     This code transform EEDL data library into 
//                     Geant4 format
//
//      Author:        V.Ivanchenko 
// 
//      Creation date: 7 November 2001
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "CLHEP/Matrix/SymMatrix.h"
#include "G4DataVector.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LinInterpolation.hh"
#include "G4LogLogInterpolation.hh"
#include "G4LinLogInterpolation.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include <string>
#include <iostream.h>
#include <stdlib.h>
#include <strstream.h>
//#include "CLHEP/Hist/HBookFile.h"
//#include "CLHEP/Hist/Histogram.h"


int main(int argc,char** argv)
{

  // -------------------------------------------------------------------
  // Setup
  G4int verbose = 0;
  /*
  HepHistogram* h1[100];
  HepHistogram* h2[100];
  hbookManager = new HBookFile(hFile, 58);
  h[0] = hbookManager->histogram("Kinetic Energy (MeV)", 50,0.,gEnergy/MeV);
  */
  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
    cout << "Input file is not specified! Exit" << endl;
    exit(1);
  }

  ifstream* fin = new ifstream();
  string fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    cout << "Input file <" << fname << "> does not exist! Exit" << endl;
    exit(1);
  }

  ofstream* fout_a = new ofstream();
  string fname1 = "eedl/br-sp.dat";
  fout_a->open(fname1.c_str(), std::ios::out|std::ios::trunc);

  ofstream* fout_b = new ofstream();
  const string fname2 = "eedl/ion-sp-";
  const string fname4 = "eedl/ion-ex-sig.dat";
  const string fname5 = "eedl/ion-ex-av.dat";

  ofstream* fout_c = new ofstream();
  fout_c->open(fname4.c_str(), std::ios::out|std::ios::trunc);
  ofstream* fout_d = new ofstream();
  fout_d->open(fname5.c_str(), std::ios::out|std::ios::trunc);

  // flags

  G4bool end     = false;
  G4bool first   = true;
  G4bool secon   = false;
  G4bool infor   = false;
  G4bool bremspec= false;
  G4bool ionis   = false;
  G4bool exsig   = false;
  G4bool exav    = false;

  //there can't be lines longer than nmax characters
  const size_t nmax = 73;
  char line[nmax]; 
  std::string line1, line2;


  // parameters

  G4double ein = 0.0;
  G4double eout = 0.0;
  //  G4double pb,fac,x,s0,sx,sx2,sx3,sx4,sx5,sx6,sy,syx,syx2,syx3,z;
  //  G4double xx,w,sz2,szx,szx2,szx3,syz;
  G4double x, z, xx, pb;
  G4double einold = 0.0;
  G4double eioold = 0.0;
  G4double eesold = 0.0;
  G4double eeaold = 0.0;
  G4DataVector e;
  G4DataVector p;
  G4DataVector dp;
  e.clear();
  p.clear();
  dp.clear();  
  G4DataVector xion;
  G4DataVector pion;
  G4DataVector tion;
  xion.clear();
  pion.clear();
  tion.clear();
  G4int nibest = 0;
  G4int nigood = 0;
  G4int nibad  = 0;
  G4double zbest = 0.01;
  G4double zgood = 0.10;

  size_t counter = 0;
  G4int Zold = -1;
  G4int shell = 0;
  G4int Z = 0; 
  G4int finalp, reactionId, dataId;
  G4double delmax = 0.0;
  G4double deltav = 0.0;
  G4double nz = 0.0;
  G4double ymax = 1.0;
  G4double zdelmax = 0.0;
  G4double dsumav = 0.0;
  G4double dsumax = 0.0;
  G4double zdeltav = 0.0;
  G4double x4max = 0.0;
  G4double znz = 0.0;
  G4double be = 0.0;

  G4double xdelmax = 0.0;
  G4double zbad = 0;
   
  // main loop 

  do {

    // read next line

    counter++;
    for(size_t ii = 0; ii < nmax; ii++) {line[ii] = ' ';}
    fin->getline( line, nmax);
    line1 = std::string("");
    line1 = std::string(line, nmax);
    if(1 < verbose) {
      cout << "Next line # " << counter << ": " << line1 << endl;
    }  

    // analize line contence

    // end of file
    if(fin->eof()) {
      end = true;
      (*fout_a) << "-2 -2" << endl;
      (*fout_b) << "-2 -2" << endl;
      (*fout_c) << "-2 -2" << endl;
      (*fout_d) << "-2 -2" << endl;

      // end of data set
    } else if(line[71] == '1') {

      first = true;
      infor = false;
      ein = 0.0;

      // first line in a data set
    } else if(first) {

      line2 = std::string(&line[0], 3);
      Z = atoi(line2.c_str());
      line2 = std::string(&line[10], 2);
      finalp = atoi(line2.c_str());
      if(0 < verbose) {
        cout << "New data for Z= " << Z
             << " finalParticleId= " << finalp 
             << endl;
      }  
      first = false;
      secon = true;

      // second line in a data set
    } else if(secon) {

      line2 = std::string(&line[0], 2);
      reactionId = atoi(line2.c_str());
      line2 = std::string(&line[2], 3);
      dataId = atoi(line2.c_str());
      if(-1 < verbose) {
        cout << "   reactionId= " << reactionId
             << "   dataId= " << dataId 
             << endl;
      }  

      // Bremsstrahlung data

      if(reactionId == 82 && (dataId == 0 || dataId == 21)) {
        verbose = 1;
        if(dataId == 21) bremspec= true;
        einold = 0.0;
        ein    = 0.0;

	// Delta-electron data

      } else if(reactionId == 81 && (dataId == 0 || dataId == 21)) {
        verbose = 1;
        eioold = 0.0;
        ein    = 0.0;
        if(dataId == 21) {
          ionis= true;
          if(Z > Zold) {
            if(Z > 1) {
              (*fout_b) << "-1 -1" << endl;
              (*fout_b) << "-2 -2" << endl;
              fout_b->close();
	    }
            shell = 0;
            Zold = Z;
            char nameChar[20] = {""};
            G4std::ostrstream ost(nameChar, 20, G4std::ios::out);
            ost << fname2 << Z << ".dat";
            string fname3(nameChar);  
            fout_b->open(fname3.c_str(), std::ios::out|std::ios::trunc);
          } else {
            shell++;
            (*fout_b) << "-1 -1" << endl;
	  }      
          if(Z <= 99) {
            be =((G4AtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy());
            if(0 < verbose) {
              G4cout << "New ioni data Z= " << Z 
                     << " shell= " << shell 
                     << " be= " << be << G4endl; 
	    }
          } else {
            ionis = false;
	  }
	}
	// Excitation

      } else if(reactionId == 83 && (dataId == 0 || dataId == 11)) {
        verbose = 1;
        eesold = 0.0;
        eeaold = 0.0;
        ein    = 0.0;
        if(dataId == 0) {
          exsig = true;
          (*fout_c) << "0 0" << endl;
	}
        if(dataId == 11) {
          exav = true;
          (*fout_d) << "0 0" << endl;
	}

      } else {
        verbose = 0;
      }

      secon = false;
      infor = true;

      // line with information in a data set

    } else if(infor) {

      if(bremspec || ionis) {
        line2 = std::string(&line[0], 15);
        ein = atof(line2.c_str());
        line2 = std::string(&line[15], 13);
        eout = atof(line2.c_str());
        line2 = std::string(&line[28], 13);
        pb = atof(line2.c_str());
        if(1 < verbose) {
          G4cout << ein << " " << eout << " " << pb << G4endl; 
	}
      } else if(exsig || exav) {
        line2 = std::string(&line[0], 15);
        ein = atof(line2.c_str());
        line2 = std::string(&line[15], 13);
        pb = atof(line2.c_str());
        if(1 < verbose) {
          G4cout << ein << " " << pb << G4endl; 
	}
      }
    }

    // handle with Bremsstrahlung data

    if(bremspec) {

      // end of sample for given energy
      if(einold > 1.0*eV && (first || ein != einold)) {

        size_t nn = e.size();
        size_t i;
        G4double x0 = 0.01;
        size_t imax = 15;
        G4double dx = 0.1;

 
        G4DataVector* eee = new G4DataVector();
        G4DataVector* ppp = new G4DataVector();
        for (i=0; i<nn; i++) {
          eee->push_back(e[i]);
          ppp->push_back(p[i]);
	}

        G4VDataSetAlgorithm* inpl = new G4LinInterpolation();
        G4VEMDataSet* dataSet = new G4EMDataSet(1, eee, ppp, inpl, 1., 1.);

        G4DataVector* eeee = new G4DataVector();
        G4DataVector* pppp = new G4DataVector();
        G4VDataSetAlgorithm* inpll = new G4LinInterpolation();
        G4double a = dataSet->FindValue(x0);
        for(i=0; i<imax; i++) {
          G4double en = dx*((G4double)i);
          if(i == 0)  en = 0.01;
          if(i == 10) en = 0.95;
          if(i == 11) en = 0.97;
          if(i == 12) en = 0.99;
          if(i == 13) en = 0.995;
          if(i == 14) en = 1.00;
          G4double pn = (dataSet->FindValue(en))/a;
          eeee->push_back(en);
          pppp->push_back(pn);
	}
        G4VEMDataSet* ndataSet = new G4EMDataSet(1, eeee, pppp, inpll, 1., 1.);
	
    
        G4double delm = 0.0;
	G4double dy, yy;
        ymax = 1.0;
        for (i=0; i<nn; i++) {
          x = e[i];
          if(x <= x0) {
            yy = 1.0;
	  } else {
	    yy = ndataSet->FindValue(x);
	  }

          if(yy > ymax) ymax = yy;

          dy = G4std::abs(p[i] - yy*a)/p[i];

          if(1 < verbose) {
            cout << "x= " << x   
                 << " p[i]= " << p[i] 
                 << " dp= " << dy 
                 << endl;
	  }
          if(dy > delm) delm = dy;
        }
   

	if(einold < 100.*MeV) {
          if(delm > delmax) delmax = delm;
          nz += 1.0;
          deltav = (deltav*(nz - 1.0) + delm)/nz;        
	}

        (*fout_a).precision(4);

        (*fout_a) << einold/MeV << " " 
                  << ymax
                  << endl;
        for(i=0; i<imax; i++) {
          (*fout_a) << (*pppp)[i] << " "; 
	}
	
        (*fout_a)  << endl;
        if(0 < verbose) {
          cout << "!!!!! e= " << einold/MeV << " x0= " << x0 
               << " dymax= " << delm
               << " dm= " << delmax
               << " da= " << deltav
               << " ymax= " << ymax 
               << endl;
	}
        einold = 0.0;
	delete dataSet;
	delete ndataSet;
      }


      // new sample
      if(ein != einold) { 
        e.clear();
        p.clear();
        dp.clear();  
        einold = ein;
      } 

      if(einold > 0.1*eV) {
        x = eout/einold;
        e.push_back(x);
        p.push_back(pb*x);
        dp.push_back(pb);

	if(1 < verbose) {
            cout << "x= " << x   
                 << " p= " << pb << " pp= " << x*pb << endl;
	}
      }
      if(first) {
        (*fout_a) << "-1 " << Z << endl;
        bremspec = false;
      }
    }
    
    if(ionis) {

      // end of sample for given energy
      size_t nn = xion.size();
      if(first || ein != eioold && (Z == 82 || Z == 27 && Z == 1)) {
        cout << "??? ein= " << ein << "; nn= " << nn << endl;
      }
      if(eioold > 1.0*eV && nn > 4 && (first || ein != eioold)) {

        G4double gam = (eioold + be)/electron_mass_c2 + 1.;
        G4double g = (2.*gam - 1.)/(gam*gam);

	G4double x1 = (0.1*eV + be)/(eioold + be);
        if(x1 > 0.5) x1 = 0.5;
        G4double x2 = 0.5; 
    
        x = xion[nn-1];
        z = 1. - g*x + (1. - g)*x*x + x*x*(1./(1. - x) - g)/(1. - x);

        G4double f1 = tion[nn-1]/z;
        G4double f2 = 0.0;
        G4double f3 = 0.0;

	// aria interpolation
        verbose = 1;
	size_t jmax1 = 3;
	size_t jmax2 = 16;
	size_t jmax  = jmax1 + jmax2 + 1;
        G4int zpr = -54;
        G4double range = (G4double)jmax;
        G4double xmax = xion[0];
        G4double tmax = tion[0];
          G4double y  = 0.0;

        for(size_t i=0; i<nn; i++) {
          tion[i] /= f1;
          if(tion[i] > tmax && i < nn-1) {
            tmax = tion[i];
            xmax = xion[i];
	  }
	} 

        G4double x3 = xmax * 10.0;
        if(x3 > 0.5)  x3 = 0.5;
        G4double dx = (x3-x1)/range;
        G4double x4 = x3 + dx;

        for(i=0; i<nn; i++) {
          x = xion[i];
          if(x < 0.499 && x >= x3 ) {
            z = 1. - g*x + (1. - g)*x*x + x*x*(1./(1. - x) - g)/(1. - x);
            f2 += (tion[i] - z)/(0.5 - x);
            f3 += 1.0/x;
	  }
	}
        if(f3 > 0.0) f2 /= f3;

        G4double zdelm = 0.0;
        G4double zdelm0= 0.0;
        G4double xdelmax0 = 0.0;
        G4double xdelmax1 = x3;
        G4bool im  = true;
        G4bool out  = false;
        size_t iiiimax = 2;

        for(size_t iiii=0; iiii<iiiimax; iiii++) {

          G4double range1 = (G4double)jmax1;
          G4double range2 = (G4double)jmax2;
          jmax  = jmax1 + jmax2 + 1;
    	  // lin
	  G4double dx1 = (xmax - x1)/range1;

	  // log
          G4double dx2 = log(x3/xmax)/range2;
          x4 = x3*exp(dx2);
	
          G4DataVector* eee = new G4DataVector();
          G4DataVector* ppp = new G4DataVector();
          G4VDataSetAlgorithm* inpl = new G4LinLogInterpolation();

          eee->push_back(x1);
          ppp->push_back(tion[0]);

          for(i=0; i<nn; i++) {
            if(xion[i] < x4) {
              eee->push_back(xion[i]);
              ppp->push_back(tion[i]);
	    }
	  } 

          G4VEMDataSet* dataSet = new G4EMDataSet(1,eee,ppp,inpl,1.,1.);

          G4DataVector* eeee = new G4DataVector();
          G4DataVector* pppp = new G4DataVector();
          G4VDataSetAlgorithm* inpll = new G4LinInterpolation();

          for(i=0; i<jmax; i++) {
	 
            if(i <= jmax1) {
	      x = x1 + dx1 * (G4double)i;
 
            } else {
          
	      x = xmax*exp(dx2 * (G4double)(i - jmax1));
	    }

            eeee->push_back(x);
            if(i < jmax-1) {
              y = dataSet->FindValue(x);
            } else {
              z = 1. - g*x + (1. - g)*x*x + x*x*(1./(1. - x) - g)/(1. - x);
              y = z + f2*(0.5 - x)/x;
	    }
	    pppp->push_back(y);
	  }
	
          G4VEMDataSet* ndataSet = new G4EMDataSet(1,eeee,pppp,inpll,1.,1.);
	    

 	  G4double dy = 0.0;
          zdelm = 0.0;
 
          x = x2;
          z = 1. - g*x + (1. - g)*x*x + x*x*(1./(1. - x) - g)/(1. - x);
          G4double zmax2 = z;
          G4double xmax2 = x2;
          G4double zmax = 0.0;
          G4bool poin = true;
          G4double yy, yyy;
          G4double dsum = 0.0;
          G4double sum = 0.0;
          G4double x0 = 0.0;
          G4double y0 = 0.0;
          G4double yy0= 0.0;

          for (i=0; i<nn; i++) {
            x = xion[nn - i - 1];
            y = tion[nn - i - 1];
            z = 1. - g*x + (1. - g)*x*x + x*x*(1./(1. - x) - g)/(1. - x);
	    yyy = z + f2*(0.5 - x)/x; 

            if(x <= x3) {
              yy = ndataSet->FindValue(x);

            } else {
              yy = yyy;
	    }
	  

            if(poin && yy > zmax2) {
              poin = false;
              xmax2 = xion[i+1];
            } else if(!poin && yy*x > zmax && x >= xmax) {
              zmax = yy*x;
	    }

            dy = 1.0 - yy/y;

            if((1 < verbose || (Z == zpr)) && iiii == 1 ) {
              cout << "x= " << x   
                   << " tion[i]= " << y 
                   << " delta= " << dy  
                   << " x1= " << x1
                   << " x2= " << xmax
                   << " x3= " << x3
                   << endl;
	    }
            dy = abs(dy);
            if(dy > zdelm) {
              zdelm = dy;
              xdelmax0 = x;
	    }
            if(im && dy > zgood) {
              im = false;
              xdelmax1 = x;
	    }
            if(i > 0) {
              G4double xsq = x*x;
              G4double xsq0= x0*x0;
              sum += (y/xsq + y0/xsq0)*(x0 - x)*0.5;
              dsum += abs((y/xsq + y0/xsq0) - (yy/xsq + yy0/xsq0))*(x0 - x)*0.5;
	    }
            yy0 = yy;
            y0  = y;
            x0  = x;
          }
   	
	
	  if(iiii == iiiimax-1 || zdelm <= zgood ) {

            out = true;	
            //if(Z == zpr) ndataSet->PrintData();
            G4String qu;
            if(dsum/sum <= zbest) {
              nibest++;
              qu = " ##BEST";
	    } else if(dsum/sum <= zgood) {
              nigood++;
              qu = " ##GOOD";
	    } else {
              nibad++;
              qu = " ##BAD";
  	    }
            if(zdelm > zdelmax) {
              zdelmax = zdelm;
              xdelmax = xdelmax0;
              zbad = Z;
	    }

            znz += 1.0;
            zdeltav = (zdeltav*(znz - 1.0) + zdelm)/znz;        

            dsum /= sum;
            if(dsum > dsumax) dsumax = dsum;
            dsumav = (dsumav*(znz - 1.0) + dsum)/znz;        

            yy0 = be + eioold;
            f2   *= yy0;
            x1   *= yy0;
            xmax *= yy0;
            x3   *= yy0;

            (*fout_b).precision(5);

            (*fout_b) << eioold/MeV << " " 
                      << f2
                      << endl;

            (*fout_b) << x1 << " ";
            (*fout_b) << xmax << " ";
            (*fout_b) << x3 << " ";


            (*fout_b).precision(4);
	
            for (i=0; i<pppp->size(); i++) {
              (*fout_b) << (*pppp)[i] << " ";
	    }
	
	
            (*fout_b) << endl;
       
            if(0 < verbose || (Z == zpr)) {
              cout << "!! e= " << eioold/MeV  
                   << " dMax= " << zdelm
                   << " dMaxAll= " << zdelmax
                   << " xMaxAll= " << xdelmax
                   << " dAv= " << zdeltav
                   << " A= " << f2 
                   << " x1= " << x1 
                   << " x2= " << xmax  
                   << " x3= " << x3 
                   << endl;
              cout << "   dsum= " << dsum
                   << " dSumMax= " << dsumax
                   << " dSumAv= " << dsumav
                   << " zmax= " << zmax 
	           << " zmax2= " << zmax2
                   << " xmax2= " << xmax2
		   << " jmax= " << jmax 
                   << " " << qu
                   << endl;
	  
              for (i=0; i<pppp->size(); i++) {
                cout << (*pppp)[i] << " ";
	      }
              cout << "  Zbad= " << zbad << endl;
	  
  	    }
	  } else if(iiii == 0) {
            zdelm0 = zdelm;
	    /*
  	    if(xdelmax1 > x3) {
              //jmax3 = 12;
              x3 = xdelmax1;
	    }
	    */
	  }
          delete dataSet;
	  delete ndataSet;
	  if(out) break;

	} 
        einold = 0.0;
      }


      // new sample
      if(ein != eioold) { 
        xion.clear();
        pion.clear();
        tion.clear();
        eioold = ein;
      } 

      if(eioold > 0.1*eV) {
        x = (eout + be)/(eioold + be);
        xion.push_back(x);
        pion.push_back(pb);
        xx = pb*x*x;
        tion.push_back(xx);
	if(1 < verbose) {
            cout << "x= " << x   
                 << " p= " << pb 
                 << " px2= " << xx
                 << " b= " << be 
                 << " ein= " << ein << endl;
	}
      }
      if(first) {
        ionis = false;
      }
    }
    if(exsig) {

      (*fout_c).precision(4);
      if(ein != eesold) eesold = ein;  
      if(eesold >= 1.0*eV && eesold <= 100.*MeV) {
        (*fout_c) << eesold/MeV << " " 
                  << pb*1.e-6
                  << endl;
        if(1 < verbose) {
          cout << "!!!!! e= " << eesold/MeV << " sig= " << pb
               << endl;
	}
      }
      if(first) {
        (*fout_c) << "-1 " << Z <<endl;
        exsig = false;
      }
    }
    if(exav) {
      (*fout_d).precision(4);
      if(ein != eeaold) eeaold = ein;  
      if(eeaold > 1.0*eV && eeaold <= 100.*MeV) {
        (*fout_d) << eeaold/MeV << " " 
                  << pb*1.e+6
                  << endl;
        if(1 < verbose) {
          cout << "!!!!! e= " << eeaold/MeV << " Eav= " << pb/MeV
               << endl;
	}
      }
      if(first) {
        (*fout_d) << "-1 " << Z << endl;
        exav = false;
      }
    }

  } while(!end && counter < 1000000);

  G4cout << "nibest = " << nibest << G4endl;
  G4cout << "nigood = " << nigood << G4endl;
  G4cout << "nibad  = " << nibad << G4endl;
  G4cout << "dMaxAll= " << zdelmax << G4endl;
  G4cout << "xMaxAll= " << xdelmax << G4endl;
  G4cout << "dAv    = " << zdeltav << G4endl;
  G4cout << "dSumMax= " << dsumax << G4endl;
  G4cout << "dSumAv = " << dsumav << G4endl;
  G4cout << "x4max = " << x4max << G4endl;
  G4cout << "Zbad   = " << zbad << G4endl;

  G4cout << "###### End of test #####" << G4endl;
}

