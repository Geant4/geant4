// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RandGlobalTest.cc,v 1.3 1999-11-23 15:00:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include "Randomize.hh"
#include "G4ios.hh"

RandEngine theRandEngine;
DRand48Engine theDRand48Engine;
RanluxEngine theRanluxEngine(19780503,4);
RanecuEngine theRanecuEngine;

//========================  CHI-SQUARED TEST  ============================

HepInt box[10], chi[8], more=1;
HepInt   prob[8] = { 1, 5, 25, 50, 75, 95, 99, 100 };
HepFloat perc[8] = { 2.09, 3.33, 5.90, 8.34, 11.39, 16.92, 21.67, 100.0 };

void even_test(HepInt seq)
{
   HepInt i,j;
   char Pause;

   for (i=0; i<8; ++i)
     chi[i]=0;

   G4cout << "\t\t\t --- CHI-SQUARED TEST ---" << G4endl;
   HepInt fac = HepInt(seq/100);

   G4cout << G4endl;

   for (HepInt k=0; k<seq; ++k) {
     for (i=0; i<10; ++i)
       box[i]=0;
     for (i=0; i<1000; ++i) {
       HepInt h = HepInt(100*G4UniformRand());
       for (j=0; j<10; ++j)
         if ((10*j <= h) && (h < 10*(j+1))) ++box[j];
     }
     HepFloat X=0;
     for (i=0; i<10; ++i) {
       HepFloat t=HepFloat(box[i]-100);
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
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   if (more == 1)
     if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;
}  // end even_test()

//======================  SERIAL CORRELATION TEST  =======================

void corr_test(HepInt seq)
{
   HepInt i,j;
   char Pause;
   HepDouble U;
   HepDouble UU;
   HepDouble UV;
   HepDouble corr = 0.;
   HepDouble stde = 0.;
   HepDouble mv   = -1./999;
   HepDouble std  = (-mv)*sqrt(HepDouble(1000*(997)/(1001)));
   HepDouble next, uni, w, mve, temp[10000];

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
   stde = sqrt((1./(seq-1))*stde);

   G4cout << "Mean value predicted : " << mv << "\t"
        << "Standard deviation predicted: " << std << G4endl
        << "Mean value effective : " << mve << "\t"
        << "Standard deviation effective: " << stde << G4endl << G4endl;
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = G4cin.get()) != '\n') exit(0);
   G4cout << G4endl;

}   // end corr_test

//===========================  SPATIAL TEST  =============================

#define PI 3.1415926

static unsigned long twotoj[16] = { 0x1L,0x2L,0x4L,0x8L,0x10L,0x20L,0x40L,
                                    0x80L,0x100L,0x200L,0x400L,0x800L,
                                    0x1000L,0x2000L,0x4000L,0x8000L };

HepFloat fnc(HepFloat x1, HepFloat x2, HepFloat x3, HepFloat x4)
{
   return  (HepFloat) sqrt(x1*x1+x2*x2+x3*x3+x4*x4);
}

void spat_test ()
{
   char Pause;
   long iy[4],jpower,k;
   HepInt i,j;
   HepFloat x1,x2,x3,x4,yprob[4];

   G4cout << "\t\t\t   --- SPATIAL TEST ---" << G4endl
        << "\tCalculates PI statistically using volume of unit n-sphere "
        << "n = 2,3,4"
        << G4endl << G4endl;
   for (i=1; i<=3; i++) iy[i]=0;
   G4cout << "\t\t\t# pts\tPI\t(4/3)PI\t(1/2)PI^2" << G4endl << G4endl;
   for (j=1; j<=15; j++) {
     for (k=twotoj[j-1]; k>=0; k--) {
       x1 = (HepFloat) G4UniformRand();
       x2 = (HepFloat) G4UniformRand();
       x3 = (HepFloat) G4UniformRand();
       x4 = (HepFloat) G4UniformRand();
       if (fnc(x1,x2,0.,0.) < 1.) ++iy[1];
       if (fnc(x1,x2,x3,0.) < 1.) ++iy[2];
       if (fnc(x1,x2,x3,x4) < 1.) ++iy[3];
     }
     jpower = twotoj[j];
     for (i=1; i<=3; i++)
       yprob[i] = (HepFloat) twotoj[i+1]*iy[i]/jpower;
     if ((jpower <= 32768) && (jpower != 0))
       G4cout << "\t\t\t" << jpower << "\t" << yprob[1] << "\t"
            << yprob[2] << "\t" << yprob[3] << G4endl;
   }
   G4cout << G4endl;
   G4cout << "\t\t\tactual" << "\t" << PI << "\t"
        << 4.*PI/3. << "\t" << .5*PI*PI << G4endl << G4endl;

   if (more != 99) {
     G4cout << "                   -----  Press <ENTER> to continue  -----";
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


void layout(HepInt seq)
{
   even_test(seq);
   corr_test(seq);
   spat_test();
}   // end layout() 


void start_test()
{
   char sel;
   HepInt seq;

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
}  // end start_test()


HepInt main() {

   init();
   start_test();
   
   return 0;
}

