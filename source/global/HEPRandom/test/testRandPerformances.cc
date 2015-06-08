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
// Test for RNG algorithms for Geant4 multi-threaded
//
// Author: A. Dotti
//
//
// History: 
//
// 5 March 2013 - First implementation
// ----------------------------------------------------------------------

#include "Randomize.hh"
#include "G4Timer.hh"
#include <iostream>
#include <map>
#include <vector>

///#include "G4MKL_engine.hh"

#define SUCCESS 1
#define FAIL 0
//Just to give some time, on my MacBook Pro Intel i7 2.6GHz: 2 min for 3M TRIALS, set a reasonable number.
//Scales linearly
//#define TRIALS 10
#define TRIALS 100000000

int engineType = 0;
//per-thread engine
static G4ThreadLocal CLHEP::HepRandomEngine *pdefaultEngine = 0 ;

void InitMe( int engineSwitch = 0 )
{
  if ( ! pdefaultEngine )
  {
      switch (engineSwitch) {
          case 0:
              pdefaultEngine = new CLHEP::Ranlux64Engine();
              break;
          case 1:
              pdefaultEngine = new CLHEP::RanecuEngine();
              break;
          case 2:
              pdefaultEngine = new CLHEP::RanluxEngine();
              break;
          case 3:
              pdefaultEngine = new CLHEP::MTwistEngine();
              break;
          case 4:
              pdefaultEngine = new CLHEP::RanshiEngine();
              break;
          //case 5:
	  //    pdefaultEngine = new G4MKL_engine();
 	  //    break;
          default:
              std::cerr<<"Non valid selection for RandomEngine, valid values from 0 (default) to 4"<<std::endl;
              abort();
              break;
      }
  } 
  G4Random::setTheEngine( pdefaultEngine );
  //G4Random::setTheSeed( 1220515164 );
}

typedef int Ret_t;
typedef Ret_t (*Test_t)(void);

//===================================
//Tests
//===================================

Ret_t testRNG()
{
  InitMe( engineType );
  CLHEP::HepRandomEngine* rng = pdefaultEngine;
  //Seed with a thread specific seed
  G4Random::setTheSeed( 123456789 );
  for ( int i = 0 ; i < TRIALS ; ++i )
    rng->flat();
  return SUCCESS;
}

G4double array[TRIALS/10];
Ret_t testRNG2()
{
  InitMe( engineType );
  CLHEP::HepRandomEngine* rng = pdefaultEngine;
  //Seed with a thread specific seed
  G4Random::setTheSeed( 123456789 );
  for ( int i = 0 ; i < TRIALS/2 ; ++i )
    rng->flatArray(2,array);
  return SUCCESS;
}

Ret_t testRNG5()
{
  InitMe( engineType );
  CLHEP::HepRandomEngine* rng = pdefaultEngine;
  //Seed with a thread specific seed
  G4Random::setTheSeed( 123456789 );
  for ( int i = 0 ; i < TRIALS/5 ; ++i )
    rng->flatArray(5,array);
  return SUCCESS;
}

Ret_t testRNG10()
{
  InitMe( engineType );
  CLHEP::HepRandomEngine* rng = pdefaultEngine;
  //Seed with a thread specific seed
  G4Random::setTheSeed( 123456789 );
  for ( int i = 0 ; i < TRIALS/10 ; ++i )
    rng->flatArray(10,array);
  return SUCCESS;
}

Ret_t testRNG20()
{
  InitMe( engineType );
  CLHEP::HepRandomEngine* rng = pdefaultEngine;
  //Seed with a thread specific seed
  G4Random::setTheSeed( 123456789 );
  for ( int i = 0 ; i < TRIALS/20 ; ++i )
    rng->flatArray(20,array);
  return SUCCESS;
}

Ret_t testRNG30()
{
  InitMe( engineType );
  CLHEP::HepRandomEngine* rng = pdefaultEngine;
  //Seed with a thread specific seed
  G4Random::setTheSeed( 123456789 );
  for ( int i = 0 ; i < TRIALS/30 ; ++i )
    rng->flatArray(30,array);
  return SUCCESS;
}

Ret_t testRNGall()
{
  InitMe( engineType );
  CLHEP::HepRandomEngine* rng = pdefaultEngine;
  //Seed with a thread specific seed
  G4Random::setTheSeed( 123456789 );
  for ( int i = 0 ; i < 10 ; ++i )
    rng->flatArray(TRIALS/10,array);
  return SUCCESS;
}
//===================================
// End Tests
//===================================

struct tt {
  Test_t aTest;
  const char* name;
  const char* mess;
};

//The list of tests to be performed
tt tests[] = { 
  {&testRNG,"flat","Use flat interface."},
  {&testRNG2,"array2","Use flatArray inteface. Array size 2"},
  {&testRNG5,"array5","Use flatArray inteface. Array size 5"},
  {&testRNG10,"array10","Use flatArray inteface. Array size 10"},
  {&testRNG20,"array20","Use flatArray inteface. Array size 20"},
  {&testRNG30,"array30","Use flatArray inteface. Array size 30"},
  {&testRNGall,"arrayAll","Use flatArray interface. Array size ALL/10"},
  {NULL,NULL,NULL}
 };

Ret_t testme( const tt& t) 
{
  G4Timer timer;
  std::cout<<std::left<<"\t  Sub-Test: ";
  std::cout<<std::left<<std::setfill(' ')<<std::setw(20)<<t.name;
  timer.Start();
  Ret_t result = t.aTest();
  timer.Stop();
  if ( result == FAIL ) 
    {
      std::cout<<std::right<<" FAILED"<<std::flush;
      std::cerr<<"Sub-Test: "<<t.name<<":"<<t.mess<<" FAILED"<<std::endl;
    }
  else 
    {
      std::cout<<std::right<<" SUCCESS\n";
    }
  std::cout<<std::right<<timer<<std::endl;
  return result;
}


int main(int,char**) 
{
    bool result = SUCCESS;
    std::cout<<"Testing Geant4 RNG"<<std::endl;
    std::cout<<"Running with number histories: "<<TRIALS<<std::endl;
    std::cout<<"Tests: "<<std::endl;
    int i = 0;
    while (true) {
        if ( tests[i].aTest == NULL ) break;
        std::cout<<std::setw(17)<<std::setfill(' ')<<std::left<<tests[i].name;
        std::cout<<std::setw(3)<<" : ";
        std::cout<<std::right<<tests[i].mess<<std::endl;
        ++i;
    }
    //0-5
    const int max = 5;//6
    for ( int t = 0 ; t < max ; ++t ) //Loop on engine types
    {
        if ( pdefaultEngine )
        {
            delete pdefaultEngine;
            pdefaultEngine = NULL;
        }
        InitMe( t );
        std::cout<<"\t================================================="<<std::endl;
        std::cout<<"\tTesting Random Engine: "<<pdefaultEngine->name()<<std::endl;
        std::cout<<"\t================================================="<<std::endl;
        G4Timer someT;
        someT.Start();
        i = 0;
        while (true) {
            if ( tests[i].aTest == NULL ) break;
            result &= testme(tests[i]);
            ++i;
        }
        someT.Stop();
        std::cout<<"\tDone, user elapsed time: "<<someT.GetUserElapsed()<<" s."<<std::endl;
    }
    if ( result == SUCCESS ) return 0;
    else return 1;
}
