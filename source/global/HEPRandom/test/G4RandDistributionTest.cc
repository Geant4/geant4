// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RandDistributionTest.cc,v 1.1 1999-01-07 16:08:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include "Randomize.hh"
#include "G4ios.hh"

HepJamesRandom theJamesEngine;
RandEngine theRandEngine;
DRand48Engine theDRand48Engine;
RanluxEngine theRanluxEngine(19780503,4);
RanecuEngine theRanecuEngine;

void init()
{
   char Pause;

   G4cout << endl << endl;
   G4cout << "---------------------------- Random shooting test -----------------------------" << endl;
   G4cout << "                             --------------------                              " << endl;
   G4cout << " >>> Random Engines available <<<" << endl << endl;
   G4cout << "   > HepJamesRandom (default)" << endl;
   G4cout << "   > Rand" << endl;
   G4cout << "   > DRand48" << endl;
   G4cout << "   > Ranlux" << endl;
   G4cout << "   > Ranecu" << endl << endl;
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = cin.get()) != '\n') exit(0);
   G4cout << endl;

}  // end init()


void layout()
{
   HepFloat m=3.0;
   const HepInt size=5;
   HepDouble vect[size];

   G4cout << " Flat ]0,1[       : " << RandFlat::shoot() << endl;
   G4cout << " Flat ]0,5[       : " << RandFlat::shoot(5) << endl;
   G4cout << " Flat ]-5,3[      : " << RandFlat::shoot(-5,3) << endl;
   G4cout << " Exp (m=1)        : " << RandExponential::shoot() << endl;
   G4cout << " Exp (m=3)        : " << RandExponential::shoot(3) << endl;
   G4cout << " Gauss (m=1)      : " << RandGauss::shoot() << endl;
   G4cout << " Gauss (m=3,v=1)  : " << RandGauss::shoot(3,1) << endl;
   G4cout << " Wigner(1,0.2)    : " << RandBreitWigner::shoot(1,0.2) << endl;
   G4cout << " Wigner(1,0.2,1)  : " << RandBreitWigner::shoot(1,0.2,1) << endl;
   G4cout << " Wigner2(1,0.2)   : " << RandBreitWigner::shootM2(1,0.2) << endl;
   G4cout << " Wigner2(1,0.2,1) : " << RandBreitWigner::shootM2(1,0.2,1) << endl;
   G4cout << " IntFlat [0,99[   : " << RandFlat::shootInt(99) << endl;
   G4cout << " IntFlat [-99,37[ : " << RandFlat::shootInt(-99,37) << endl;
   G4cout << " Poisson (m=3.0)  : " << RandPoisson::shoot(m) << endl;
   G4cout << endl;
   G4cout << " Shooting an array of 5 flat numbers ..." << endl << endl;
   RandFlat::shootArray(size,vect);
   for ( HepInt i=0; i<size; ++i )
     G4cout << " " << vect[i];
   G4cout << endl << endl;
}   // end layout() 

void dist_layout()
{
   HepFloat m=3.0;
   const HepInt size=5;
   HepDouble vect[size];
   char Pause;

   HepJamesRandom aJamesEngine;
   RandEngine aRandEngine;
   DRand48Engine aDRand48Engine;
   RanluxEngine aRanluxEngine(19780503,4);
   RanecuEngine aRanecuEngine;

   RandFlat aFlatObj(aJamesEngine);
   RandExponential anExponentialObj(aRandEngine);
   RandGauss aGaussObj(aDRand48Engine);
   RandBreitWigner aBreitObj(aRanluxEngine);
   RandPoisson aPoissonObj(aRanecuEngine);

   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = cin.get()) != '\n') exit(0);
   G4cout << "-------------------- Shooting test on distribution objects --------------------" << endl;
   G4cout << endl;
   G4cout << " Flat ]0,1[       : " << aFlatObj.fire() << endl;
   G4cout << " Flat ]0,5[       : " << aFlatObj.fire(5) << endl;
   G4cout << " Flat ]-5,3[      : " << aFlatObj.fire(-5,3) << endl;
   G4cout << " Exp (m=1)        : " << anExponentialObj.fire() << endl;
   G4cout << " Exp (m=3)        : " << anExponentialObj.fire(3) << endl;
   G4cout << " Gauss (m=1)      : " << aGaussObj.fire() << endl;
   G4cout << " Gauss (m=3,v=1)  : " << aGaussObj.fire(3,1) << endl;
   G4cout << " Wigner(1,0.2)    : " << aBreitObj.fire(1,0.2) << endl;
   G4cout << " Wigner(1,0.2,1)  : " << aBreitObj.fire(1,0.2,1) << endl;
   G4cout << " Wigner2(1,0.2)   : " << aBreitObj.fireM2(1,0.2) << endl;
   G4cout << " Wigner2(1,0.2,1) : " << aBreitObj.fireM2(1,0.2,1) << endl;
   G4cout << " IntFlat [0,99[   : " << aFlatObj.fireInt(99) << endl;
   G4cout << " IntFlat [-99,37[ : " << aFlatObj.fireInt(-99,37) << endl;
   G4cout << " Poisson (m=3.0)  : " << aPoissonObj.fire(m) << endl;
   G4cout << endl;
   G4cout << " Shooting an array of 5 flat numbers ..." << endl << endl;
   aFlatObj.fireArray(size,vect);
   for ( HepInt i=0; i<size; ++i )
     G4cout << " " << vect[i];
   G4cout << endl << endl;
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = cin.get()) != '\n') exit(0);
}   // end dist_layout() 

void user_layout()
{
   HepFloat m=3.0;
   const HepInt size=5;
   HepDouble vect[size];
   char sel;
   HepRandomEngine* anEngine;

   G4cout << endl << endl;
   G4cout << "-------------------- Shooting test skeeping the generator ---------------------" << endl;
   G4cout << endl;
   G4cout << " >>> Select a Random Engine <<<" << endl << endl;
   G4cout << "   1. HepJamesRandom (default)" << endl;
   G4cout << "   2. Rand" << endl;
   G4cout << "   3. DRand48" << endl;
   G4cout << "   4. Ranlux" << endl;
   G4cout << "   5. Ranecu" << endl << endl;
   G4cout << " > ";
   cin >> sel;
   while ((sel!='1')&&(sel!='2')&&(sel!='3')&&(sel!='4')&&(sel!='5')) {
     G4cout << endl << " >>> Choice not legal !!  [1..5]<<<" << endl;
     cin >> sel;
   }

   switch (sel) {
     case '1':
       anEngine = &theJamesEngine;
       break;
     case '2':
       anEngine = &theRandEngine;
       break;
     case '3':
       anEngine = &theDRand48Engine;
       break;
     case '4':
       anEngine = &theRanluxEngine;
       break;
     case '5':
       anEngine = &theRanecuEngine;
       break;

     default:
       anEngine = &theJamesEngine;
       break;
   }
   G4cout << endl;
 
   G4cout << " Flat ]0,1[       : " << RandFlat::shoot(anEngine) << endl;
   G4cout << " Flat ]0,5[       : " << RandFlat::shoot(anEngine,5) << endl;
   G4cout << " Flat ]-5,3[      : " << RandFlat::shoot(anEngine,-5,3) << endl;
   G4cout << " Exp (m=1)        : " << RandExponential::shoot(anEngine) << endl;
   G4cout << " Exp (m=3)        : " << RandExponential::shoot(anEngine,3) << endl;
   G4cout << " Gauss (m=1)      : " << RandGauss::shoot(anEngine) << endl;
   G4cout << " Gauss (m=3,v=1)  : " << RandGauss::shoot(anEngine,3,1) << endl;
   G4cout << " Wigner(1,0.2)    : " << RandBreitWigner::shoot(anEngine,1,0.2) << endl;
   G4cout << " Wigner(1,0.2,1)  : " << RandBreitWigner::shoot(anEngine,1,0.2,1) << endl;
   G4cout << " Wigner2(1,0.2)   : " << RandBreitWigner::shootM2(anEngine,1,0.2) << endl;
   G4cout << " Wigner2(1,0.2,1) : " << RandBreitWigner::shootM2(anEngine,1,0.2,1) << endl;
   G4cout << " IntFlat [0,99[   : " << RandFlat::shootInt(anEngine,99) << endl;
   G4cout << " IntFlat [-99,37[ : " << RandFlat::shootInt(anEngine,-99,37) << endl;
   G4cout << " Poisson (m=3.0)  : " << RandPoisson::shoot(anEngine,m) << endl;
   G4cout << endl;
   G4cout << " Shooting an array of 5 flat numbers ..." << endl << endl;
   RandFlat::shootArray(anEngine,size,vect);
   for ( HepInt i=0; i<size; ++i )
     G4cout << " " << vect[i];
   G4cout << endl << endl;
}   // end layout() 

void start_test()
{
   char Pause;

   G4cout << "-------------------------  Test on HepJamesRandom  ----------------------------" << endl;
   G4cout << endl;
   layout();
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = cin.get()) != '\n') exit(0);
   G4cout << endl;
   G4cout << "---------------------------  Test on RandEngine  ------------------------------" << endl;
   G4cout << endl;
   HepRandom::setTheEngine(&theRandEngine);
   layout();
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = cin.get()) != '\n') exit(0);
   G4cout << endl;
   G4cout << "-------------------------  Test on DRand48Engine  -----------------------------" << endl;
   G4cout << endl;
   HepRandom::setTheEngine(&theDRand48Engine);
   layout();
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = cin.get()) != '\n') exit(0);
   G4cout << endl;
   G4cout << "---------------------  Test on RanluxEngine (luxury 4) ------------------------" << endl;
   G4cout << endl;
   HepRandom::setTheEngine(&theRanluxEngine);
   layout();
   G4cout << "                   -----  Press <ENTER> to continue  -----";
   if ( (Pause = cin.get()) != '\n') exit(0);
   G4cout << endl;
   G4cout << "--------------------------  Test on RanecuEngine ------------------------------" << endl;
   G4cout << endl;
   HepRandom::setTheEngine(&theRanecuEngine);
   layout();
   dist_layout();
   user_layout();
}  // end start_test()


HepInt main() {

   init();
   start_test();
   
   return 0;
}

