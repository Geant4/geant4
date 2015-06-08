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
//
// $Id$
//
// 
// ----------------------------------------------------------------------
#include "Randomize.hh"
#include "G4ios.hh"
using namespace CLHEP;
//RandEngine theRandEngine;
//DRand48Engine theDRand48Engine;
RanluxEngine theRanluxEngine(19780503,4);
RanecuEngine theRanecuEngine;
#include "G4MKL_engine.hh"
G4MKL_engine theMKLengine;

//========================  CHI-SQUARED TEST  ============================

int box[10], chi[8], more=1;
int   prob[8] = { 1, 5, 25, 50, 75, 95, 99, 100 };
float perc[8] = { 2.09, 3.33, 5.90, 8.34, 11.39, 16.92, 21.67, 100.0 };

void even_test(int seq)
{
   int i,j;
   char Pause;

   for (i=0; i<8; ++i)
     chi[i]=0;

   G4cout << "\t\t\t --- CHI-SQUARED TEST ---" << G4endl;
   int fac = int(seq/100);

   G4cout << G4endl;

   for (int k=0; k<seq; ++k) {
     for (i=0; i<10; ++i)
       box[i]=0;
     for (i=0; i<1000; ++i) {
       int h = int(100*G4UniformRand());
       for (j=0; j<10; ++j)
         if ((10*j <= h) && (h < 10*(j+1))) ++box[j];
     }
     float X=0;
     for (i=0; i<10; ++i) {
       float t=float(box[i]-100);
       X += t*t;
     }
     X /= 100;
     if (k == (seq%100)) {
       G4cout << "\tDistribution of a single run (rand_valuex100) ..."
            << G4endl << G4endl;
       G4cout << "\t\t 1   10  20  30  40  50  60  70  80  90  100" << G4endl;
       G4cout << "\t\t";
       for (j=0; j<10; ++j)
         G4cout << box[j] << "  ";
       G4cout << G4endl;
       G4cout << "\t\t\t     Chi-squared = " << X << G4endl << G4endl;
     }
     for (i=0; i<8; ++i)
       if (X <= perc[i]) ++chi[i];
   }
   G4cout << "\tDistribution of Chi-squared for " << seq << " sequences ..."
        << G4endl << G4endl;
   G4cout << "\t\t\t       %" << "\tChi-sq " << "\tExpected" << G4endl;
   for (j=0; j<8; ++j)
     G4cout << "\t\t\tX <= " << perc[j] << ":\t" << chi[j]/fac
          << "\t" << prob[j] << G4endl;
   G4cout << G4endl;
   G4cout << G4endl;
   G4cout << "                   -----  Press <ENTER> to continue  -----"<<G4endl;
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   if (more == 1)
     if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;
}  // end even_test()

//======================  SERIAL CORRELATION TEST  =======================

void corr_test(int seq)
{
   int i,j;
   char Pause;
   double U;
   double UU;
   double UV;
   double corr = 0.;
   double stde = 0.;
   double mv   = -1./999;
   double std  = (-mv)*std::sqrt(double(1000*(997)/(1001)));
   double next, uni, w, mve, temp[10000];

   G4cout << "\t\t\t--- SERIAL CORRELATION TEST ---" << G4endl << G4endl;

   for (j=0; j<seq; ++j) {
     U = 0.; UU = 0.; UV = 0.;
     uni = G4UniformRand();
     w = uni;
     for (i=1; i<1000; ++i) {
       next = G4UniformRand();
       U  += next;
       UU += next*next;
       UV += uni*next;
       uni = next;
     }
     U  += w;
     UU += w*w;
     UV += uni*w;
     temp[j] = (1000*UV-U*U)/(1000*UU-U*U);
     corr += temp[j];
   }

   mve  = corr/seq;
   for (j=0; j<seq; ++j) {
     temp[j] -= mve;
     stde += temp[j]*temp[j];
   }
   stde = std::sqrt((1./(seq-1))*stde);

   G4cout << "Mean value predicted : " << mv << "\t"
        << "Standard deviation predicted: " << std << G4endl
        << "Mean value effective : " << mve << "\t"
        << "Standard deviation effective: " << stde << G4endl << G4endl;
   G4cout << "                   -----  Press <ENTER> to continue  -----"<<G4endl;
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;

}   // end corr_test

//===========================  SPATIAL TEST  =============================

#define PI 3.1415926

static unsigned long twotoj[16] = { 0x1L,0x2L,0x4L,0x8L,0x10L,0x20L,0x40L,
                                    0x80L,0x100L,0x200L,0x400L,0x800L,
                                    0x1000L,0x2000L,0x4000L,0x8000L };

float fnc(float x1, float x2, float x3, float x4)
{
   return  (float) std::sqrt(x1*x1+x2*x2+x3*x3+x4*x4);
}

void spat_test ()
{
   char Pause;
   long iy[4],jpower,k;
   int i,j;
   float x1,x2,x3,x4,yprob[4];

   G4cout << "\t\t\t   --- SPATIAL TEST ---" << G4endl
        << "\tCalculates PI statistically using volume of unit n-sphere "
        << "n = 2,3,4"
        << G4endl << G4endl;
   for (i=1; i<=3; i++) iy[i]=0;
   G4cout << "\t\t\t# pts\tPI\t(4/3)PI\t(1/2)PI^2" << G4endl << G4endl;
   for (j=1; j<=15; j++) {
     for (k=twotoj[j-1]; k>=0; k--) {
       x1 = (float) G4UniformRand();
       x2 = (float) G4UniformRand();
       x3 = (float) G4UniformRand();
       x4 = (float) G4UniformRand();
       if (fnc(x1,x2,0.,0.) < 1.) ++iy[1];
       if (fnc(x1,x2,x3,0.) < 1.) ++iy[2];
       if (fnc(x1,x2,x3,x4) < 1.) ++iy[3];
     }
     jpower = twotoj[j];
     for (i=1; i<=3; i++)
       yprob[i] = (float) twotoj[i+1]*iy[i]/jpower;
     if ((jpower <= 32768) && (jpower != 0))
       G4cout << "\t\t\t" << jpower << "\t" << yprob[1] << "\t"
            << yprob[2] << "\t" << yprob[3] << G4endl;
   }
   G4cout << G4endl;
   G4cout << "\t\t\tactual" << "\t" << PI << "\t"
        << 4.*PI/3. << "\t" << .5*PI*PI << G4endl << G4endl;

   if (more != 99) {
     G4cout << "                   -----  Press <ENTER> to continue  -----"<<G4endl;
     if ( (Pause = G4cin.get()) != '\n') exit(0);
     G4cout << G4endl;
   }
}  // end spat_test()

//=============================  GLOBAL  =================================

void init()
{
   G4cout << G4endl << G4endl;
   G4cout << "-------------------------- Random distribution test ---------------------------" << G4endl;
   G4cout << "                           ------------------------                            " << G4endl;
   G4cout << " >>> Random Engines available <<< " << G4endl << G4endl;
   G4cout << "   > HepJamesRandom (default)" << G4endl;
   G4cout << "   > Rand" << G4endl;
   G4cout << "   > DRand48" << G4endl;
   G4cout << "   > Ranlux" << G4endl;
   G4cout << "   > Ranecu" << G4endl << G4endl;
   G4cout << " >>> Tests performed <<< " << G4endl << G4endl;
   G4cout << "   > Even distribution test" << G4endl;
   G4cout << "   > Serial correlation test" << G4endl;
   G4cout << "   > Spatial test" << G4endl;
   G4cout << G4endl << G4endl;

}  // end init()


void layout(int seq)
{
   even_test(seq);
   corr_test(seq);
   spat_test();
}   // end layout() 


void start_test()
{
   char sel;
   int seq;

   G4cout << " Select the number of random values sequences:" << G4endl;
   G4cout << "\t   a - 100 sequences of 1000 random values" << G4endl;
   G4cout << "\t   b - 1000 sequences of 1000 random values" << G4endl;
   G4cout << "\t   c - 10000 sequences of 1000 random values" << G4endl;
   G4cout << "   > ";
   G4cin >> sel;
   if ((sel!='a')&&(sel!='b')&&(sel!='c'))
     exit(0);
   switch (sel) {
     case 'a':
       seq = 100;
       break;
     case 'b':
       seq = 1000;
       break;
     case 'c':
       seq = 10000;
       break;
     default:
       seq = 100;
       break;
   }

   G4cout << G4endl << G4endl;
   G4cout << "-------------------------  Test on HepJamesRandom  ----------------------------" << G4endl;
   G4cout << G4endl;
   layout(seq);
   G4cout << G4endl << G4endl;
   more = 0;

   /*
   G4cout << "---------------------------  Test on RandEngine  ------------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theRandEngine);
   layout(seq);
   G4cout << G4endl << G4endl;
   G4cout << "-------------------------  Test on DRand48Engine  -----------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theDRand48Engine);
   layout(seq);
   G4cout << G4endl << G4endl;
   */
   G4cout << "---------------------  Test on RanluxEngine (luxury 4) ------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theRanluxEngine);
   layout(seq);
   G4cout << G4endl << G4endl;
   more = 99;
   G4cout << "--------------------------  Test on RanecuEngine ------------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theRanecuEngine);
   layout(seq);

   G4cout << "--------------------------  Test on MKL G4 wrapper ------------------------------" << G4endl;
   G4cout << G4endl;
   HepRandom::setTheEngine(&theMKLengine);
   layout(seq);


}  // end start_test()



int main() {

   init();
   start_test();
   
   return 0;
}

