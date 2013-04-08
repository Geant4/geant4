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
#include "G4Threading.hh"
#include "G4Timer.hh"
#include "G4AutoLock.hh"
#include <iostream>
#include <map>
#include <vector>
#include <list>

#define SUCCESS 1
#define FAIL 0
//Just to give some time, on my MacBook Pro Intel i7 2.6GHz: 2 min for 3M TRIALS, set a reasonable number.
//Scales linearly
#define TRIALS 10000000
#define NTHREADS 4
#define STATFILE "/tmp/status.rand"


int engineType = 0 ;

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
          default:
              std::cerr<<"Non valid selection for RandomEngine, valid values from 0 (default) to 4"<<std::endl;
              abort();
              break;
      }
  }
  G4Random::setTheEngine( pdefaultEngine );
  //G4Random::setTheSeed( 1220515164 );
}

G4Mutex mutex = G4MUTEX_INITIALIZER;

//Argument to thread function
struct arg_t {
  long seed;
  long histories;
};

std::vector<std::list<uint32_t>*> outputLists;


//Simply run and store final value
G4ThreadFunReturnType genRandomNumbers(G4ThreadFunArgType arg)
{
  InitMe( engineType );
  arg_t* conf = (arg_t*)arg;
  long seed = conf->seed;
    std::list<uint32_t> * values = new std::list<uint32_t>;
  G4Random::setTheSeed( seed );
  for ( int i = 0 ; i < conf->histories ; ++i )
      values->push_back(G4RandFlat::shootInt( long(0), long(UINT32_MAX) ) );
  G4AutoLock l(&mutex);
    std::cout<<"Done thread: "<<G4GetPidId()<<std::endl;
  outputLists.push_back(values);
  return (G4ThreadFunReturnType)NULL;
}

void generate()
{
  G4Thread threads[NTHREADS];
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
        arg_t myarg = { G4RandFlat::shootInt( long(0),long(LONG_MAX) ) , TRIALS };
      G4THREADCREATE( &(threads[i]) ,  genRandomNumbers , &myarg );
    }
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      G4THREADJOIN( threads[i] );
    }
    std::cout<<"Generation done"<<std::endl;
}


#include <fstream>
#include <iomanip>
int main(int argc,char** argv)
{
    std::cout<<"Testing Geant4 RNG (works only in MT build)"<<std::endl;
    std::cout<<"Preparing output ascii file for Marsaglia battery test."<<std::endl;
    std::cout<<"Running with "<<NTHREADS<<" threads, total number histories: "<<TRIALS*NTHREADS<<std::endl;
    if ( argc>1 )
        InitMe( atoi(argv[1]) );
    else
        InitMe();
    std::cout<<"Using Engine: "<<pdefaultEngine->name()<<std::endl;
    G4Random::setTheSeed( 1234567890 );
    G4Timer someT;
    someT.Start();
    generate();
    std::cout<<"Writing out in data.bin"<<std::endl;
    std::ofstream myFile,myFile2;
    myFile2.open ("data.bin", std::ofstream::out|std::ofstream::binary);
    //myFile.open ("data.txt", std::ofstream::out); //TXT format for asc2bin program
    int i = 0;
    while ( ! outputLists.empty() )
    {
        std::list<uint32_t>* nums = outputLists.back();
        for ( std::list<uint32_t>::const_iterator it = nums->begin(); it != nums->end() ; ++it )
        {
            ++i;
            uint32_t anum = *it;
            //myFile<<std::setfill('0')<<std::setw(8)<<std::hex<<anum;
            myFile2.write(reinterpret_cast<char*>(&anum),sizeof(uint32_t));
            if ( i % 10 == 0 )
            {
                if ( i % 1000000 == 0 )
                    std::cout<<std::dec<<i<<": 0x"<<std::setfill('0')<<std::setw(8)<<std::hex<<anum<<" ("<<std::dec<<anum<<")"<<std::endl;
                myFile<<std::endl;
            }
        }
        delete nums;
        outputLists.pop_back();
    }
    //myFile.close();
    myFile2.close();
    return 0;
}
