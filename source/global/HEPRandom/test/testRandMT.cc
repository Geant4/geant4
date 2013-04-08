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

#define SUCCESS 1
#define FAIL 0
//Just to give some time, on my MacBook Pro Intel i7 2.6GHz: 2 min for 3M TRIALS, set a reasonable number.
//Scales linearly
//#define TRIALS 10
#define TRIALS 10000000
#define NTHREADS 4
#define STATFILE "/tmp/status.rand"

G4Mutex mutex = G4MUTEX_INITIALIZER;

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
          default:
              std::cerr<<"Non valid selection for RandomEngine, valid values from 0 (default) to 4"<<std::endl;
              abort();
              break;
      }
  } 
  G4Random::setTheEngine( pdefaultEngine );
  //G4Random::setTheSeed( 1220515164 );
}

std::map<G4Pid_t,double> finalRandomValue;

//RNG function, e.g. G4RandFlat::soot() typedef
typedef double (*rng_t)(void);

//======================================
// Functions executed in threads
//======================================

//Argument to thread function
struct arg_t {
  long seed;
  long histories;
  long halfway;
  rng_t rng;
};



//Simply run and store final value
G4ThreadFunReturnType mythreadfunc(G4ThreadFunArgType arg)
{
  InitMe( engineType );
  arg_t* conf = (arg_t*)arg;
  long seed = conf->seed;
  rng_t rng = conf->rng;
  double val;
  G4Random::setTheSeed( seed );
  for ( int i = 0 ; i < conf->histories ; ++i )
    val=rng();
  G4AutoLock l(&mutex);
  //std::cout<<val<<std::endl;
  finalRandomValue[ G4GetPidId() ] = val;
  return (G4ThreadFunReturnType)NULL;
}

//Weak repro per thread
G4ThreadFunReturnType mythreadfunc2(G4ThreadFunArgType arg)
{
  InitMe( engineType );
  arg_t* conf = (arg_t*)arg;
  long seed = conf->seed;
  rng_t rng = conf->rng;
  G4Random::setTheSeed( seed );
  char buf[200];
  sprintf(buf,"%s.%d",STATFILE,G4GetPidId());
  G4Random::saveEngineStatus( buf );
  double val1,val2;
  for ( int i = 0 ; i < conf->histories ; ++i )
    val1=rng();
  //Reseed and restart
  G4Random::restoreEngineStatus( buf );
  for ( int i = 0 ; i < conf->histories ; ++i )
    val2=rng();
  G4AutoLock l(&mutex);
  finalRandomValue[ G4GetPidId() ] = val1-val2;
  return (G4ThreadFunReturnType)NULL;
}

//Strong repro per thread
G4ThreadFunReturnType mythreadfunc3(G4ThreadFunArgType arg)
{
  InitMe( engineType );
  arg_t* conf = (arg_t*)arg;
  rng_t rng = conf->rng;
  G4Random::setTheSeed( conf->seed );
  char buf[200];
  std::vector<double> run1(conf->histories);
  for ( int i = 0 ; i < conf->histories ; ++i )
    {
      run1[i]=rng();
      if ( i == conf->halfway ) 
	{
	  sprintf(buf,"%s.%d",STATFILE,G4GetPidId());
	  G4Random::saveEngineStatus(buf);
	  //std::cout<<buf<<std::endl;
	}
    }
  //Reseed and restart
  G4Random::restoreEngineStatus( buf );
  //Now generate events from mid+1 to max
  //Random numbers should be the same as in run1
  for ( int i = (conf->halfway)+1 ; i < conf->histories ; ++i )
    {
      double val = rng();
      if ( val != run1[i] )
	{
	  G4AutoLock l(&mutex);
	  finalRandomValue[ G4GetPidId() ] = val-run1[i];
	  std::cerr<<"Error for run1["<<i<<"]="<<run1[i]<<" , value="<<val<<" (ThreadID:"<<G4GetPidId()<<")"<<std::endl;
	}
    } 
  return (G4ThreadFunReturnType)NULL;
}

//MT Check no memory (RNG does not depends on its own history)
G4ThreadFunReturnType mythreadfunc4(G4ThreadFunArgType arg)
{
  InitMe( engineType );
  arg_t* conf = (arg_t*)arg;
  rng_t rng = conf->rng;
  //Seed with a thread specific seed
  G4Random::setTheSeed( G4GetPidId() );
  for ( int i = 0 ; i < conf->halfway ; ++i )
      rng();
  //Now reseed with common seed
  G4Random::setTheSeed( conf->seed );
  double val;
  //Continue
  for ( int i = (conf->halfway)+1 ; i < conf->histories ; ++i )
    val = rng();
  G4AutoLock l(&mutex);
  finalRandomValue[ G4GetPidId() ] = val;
  return (G4ThreadFunReturnType)NULL;
}

//MT Check no memory (RNG does not depends on its own history)
//Similar to previous but now common part is at the beginning
G4ThreadFunReturnType mythreadfunc5(G4ThreadFunArgType arg)
{
  InitMe( engineType );
  arg_t* conf = (arg_t*)arg;
  rng_t rng = conf->rng;
  //Seed with a thread specific seed
  G4Random::setTheSeed( conf->seed );
  for ( int i = 0 ; i < conf->halfway ; ++i )
      rng();
  //Now reseed with thread specific seed
  G4Random::setTheSeed( G4GetPidId() );
  double val;
  //Continue
  for ( int i = (conf->halfway)+1 ; i < conf->histories ; ++i )
    val = rng();
  G4AutoLock l(&mutex);
  finalRandomValue[ G4GetPidId() ] = val;
  return (G4ThreadFunReturnType)NULL;
}


typedef int Ret_t;
typedef Ret_t (*Test_t)(void);

//===================================
//Tests
//===================================
Ret_t basicTest1(void) 
{
  InitMe( engineType );
  G4RandFlat::shoot();
  G4RandBit::shoot();
  G4RandGauss::shoot();
  G4RandGamma::shoot();
  G4RandExponential::shoot();
  for ( int n = 0 ; n < 100000 ; n++ )
    G4RandFlat::shoot();
  return SUCCESS;
}

//Weak Reproducibility, simple, sequential
Ret_t sequence(rng_t rng)
{
    //if ( pdefaultEngine ) delete pdefaultEngine;
    //pdefaultEngine = 0;
  InitMe( engineType );
  long seed = 1234567890;
  long max = TRIALS;
  double value1,value2;
  //Run 1
  G4Random::setTheSeed( seed );
  G4Random::saveEngineStatus( STATFILE );
  for ( long n = 0 ; n< max ; ++n )
    value1 = rng();
  //Run 2
   G4Random::restoreEngineStatus( STATFILE );
  //G4Random::setTheSeed( seed ); //This is not enought for MTwistEngine, need full status
  for ( long n = 0 ; n< max ; ++n )
    value2 = rng();
  if ( value1 == value2 )
    return SUCCESS;
  else
  {
      std::cout<<std::endl<<value1<<" "<<value2<<std::endl;
      return FAIL;
  }
}
Ret_t sequence()
{
  bool success = sequence(&G4RandFlat::shoot);
  success &= sequence(&G4RandBit::shoot);
  success &= sequence(&G4RandGamma::shoot);
  success &= sequence(&G4RandGauss::shoot);
  success &= sequence(&G4RandExponential::shoot);
  return success;
}

//Strong Reproducibility, sequential, for a given distribution
Ret_t repro(rng_t rng )
{
  InitMe( engineType );
  long max = TRIALS;
  long mid = TRIALS/2;
    std::vector<double> run1(TRIALS);
  //Run1: max random nums 
  for ( int i = 0 ; i < max ; ++i )
    {
      run1[i]=rng(); //Shoot integers, so can be used as seeds
      if ( i == mid ) 
	{
	  //Save random engine status
	  G4Random::saveEngineStatus( STATFILE );
	}
    }
  //Reset seed starting from event mid
  G4Random::restoreEngineStatus( STATFILE );
  //Now generate events from mid+1 to max
  //Random numbers should be the same as in run1
  bool success = SUCCESS;
  for ( int i = mid+1 ; i < max ; ++i )
    {
      double val = rng();
      if ( val != run1[i] )
      {
        //std::cerr<<"Error for run1["<<i<<"]="<<run1[i]<<" , value="<<val<<std::endl;
        success = FAIL;
      }
    }
  return success;
}

//Test all distributions
Ret_t repro()
{
  bool success = repro(&G4RandFlat::shoot);
  success &= repro(&G4RandBit::shoot);
  success &= repro(&G4RandGamma::shoot);
  success &= repro(&G4RandGauss::shoot);
  success &= repro(&G4RandExponential::shoot);
  return success;
}

//Simple MT test: start some threads all with same initial seed.
//Check final number is the same for all threads.
//Specific distribution
Ret_t threads1(rng_t rng)
{
  finalRandomValue.clear();
  G4Thread threads[NTHREADS];
  long seed = 123456789;
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      arg_t myarg = { seed , TRIALS , 0 , rng};
      G4THREADCREATE( &(threads[i]) ,  mythreadfunc , &myarg );
    }
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      G4THREADJOIN( threads[i] );
    }
  std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
  double val = (*it).second;
  ++it;
  for ( ; it != finalRandomValue.end() ; ++it )
    {
      if ( (*it).second != val ) return FAIL;
    }
  return SUCCESS;
}

//All distributions
Ret_t threads1()
{
  bool success = threads1( &G4RandFlat::shoot );
  success &= threads1(&G4RandBit::shoot);
  success &= threads1(&G4RandGamma::shoot);
  success &= threads1(&G4RandGauss::shoot);
  success &= threads1(&G4RandExponential::shoot);
  return success;
}

//MT Weak Reproducibility
Ret_t threads2(rng_t rng)
{
  finalRandomValue.clear();
  G4Thread threads[NTHREADS];
  InitMe( engineType );
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      arg_t myarg = { G4RandFlat::shootInt(100000) , TRIALS , 0 , rng};
      G4THREADCREATE( &(threads[i]) ,  mythreadfunc2 , &myarg );
      
    }
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      G4THREADJOIN( threads[i] );
    }
  std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
  for ( ; it != finalRandomValue.end() ; ++it )
    {
      if ( (*it).second != 0 ) return FAIL;
    }
  return SUCCESS;
}
//All distributions
Ret_t threads2()
{
  bool success = threads2( &G4RandFlat::shoot );
  success &= threads2(&G4RandBit::shoot);
  success &= threads2(&G4RandGamma::shoot);
  success &= threads2(&G4RandGauss::shoot);
  success &= threads2(&G4RandExponential::shoot);
  return success;
}


//MT Strong reproducibility
Ret_t threads3(rng_t rng)
{
  finalRandomValue.clear();
  G4Thread threads[NTHREADS];
  InitMe( engineType );
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      arg_t myarg = { G4RandFlat::shootInt(100000) , TRIALS , TRIALS/2, rng};
      G4THREADCREATE( &(threads[i]) ,  mythreadfunc3 , &myarg );
    }
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      G4THREADJOIN( threads[i] );
    }
  std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
  for ( ; it != finalRandomValue.end() ; ++it )
    {
      if ( (*it).second != 0 ) return FAIL;
    }
  return SUCCESS;
}
//All distributions
Ret_t threads3()
{
  bool success = threads3( &G4RandFlat::shoot );
  success &= threads3(&G4RandBit::shoot);
  success &= threads3(&G4RandGamma::shoot);
  success &= threads3(&G4RandGauss::shoot);
  success &= threads3(&G4RandExponential::shoot);
  return success;
}

//Check memory: reseed halfway, final random seed does not depend
//on initial seed 
Ret_t threads4(rng_t rng)
{
  finalRandomValue.clear();
  G4Thread threads[NTHREADS];
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      arg_t myarg = { G4RandFlat::shootInt(100000) , TRIALS , TRIALS/2, rng};
      G4THREADCREATE( &(threads[i]) ,  mythreadfunc4 , &myarg );
    }
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      G4THREADJOIN( threads[i] );
    }
  std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
  double val = (*it).second;
  ++it;
  for ( ; it != finalRandomValue.end() ; ++it )
    {
      if ( (*it).second != val ) return FAIL;
    }
  return SUCCESS;
}

//All distributions
Ret_t threads4()
{
  bool success = threads4( &G4RandFlat::shoot );
  success &= threads4(&G4RandBit::shoot);
  success &= threads4(&G4RandGamma::shoot);
  success &= threads4(&G4RandGauss::shoot);
  success &= threads4(&G4RandExponential::shoot);
  return success;
}

//Like previous but other way around. Halfway reseed with per-thread
//seed, final random number not found
Ret_t threads5(rng_t rng)
{
  finalRandomValue.clear();
  G4Thread threads[NTHREADS];
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      arg_t myarg = { G4RandFlat::shootInt(100000) , TRIALS , TRIALS/2, rng};
      G4THREADCREATE( &(threads[i]) ,  mythreadfunc5 , &myarg );
    }
  for ( int i = 0 ; i < NTHREADS ; ++i )
    {
      G4THREADJOIN( threads[i] );
    }
  std::map<double,bool> alreadyFound;
  std::map<G4Pid_t,double>::const_iterator it = finalRandomValue.begin();
  for ( ; it != finalRandomValue.end() ; ++it )
    {
      if ( alreadyFound.find( (*it).second ) == alreadyFound.end() ) //Not found
          alreadyFound[(*it).second]=true;
      else
          return FAIL;

    }
  return SUCCESS;
}

//All distributions
Ret_t threads5()
{
  bool success = threads5( &G4RandFlat::shoot );
  success &= threads5(&G4RandBit::shoot);
  success &= threads5(&G4RandGamma::shoot);
  success &= threads5(&G4RandGauss::shoot);
  success &= threads5(&G4RandExponential::shoot);
  return success;
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
  {&basicTest1 ,"Simple","Initialize and Generate few numbers, check no crashes."} ,
  {&sequence ,"WeakSequential","Weak Reproducibility: run twice from same initial seed. Check final values are the same."},
  {&repro ,"StrongSequential","Strong Reproducibility: run twice with reseed half way. Check final values are the same."},
  {&threads1,"MTSimple","MT simple: start threads with same random seed, check final values are the same."},
  {&threads2,"MTWeak","MT Weak Reproducibility: Same as in Weak reproducibility for sequential."},
  {&threads3,"MTStrong","MT Strong Reproducibility: same as Strong reproducibility for sequential."},
  {&threads4,"MTNoMem1","MT no-memory: threads are reseeded with common seed halfway, check final values area the same."},
  {&threads5,"MTNoMem2","MT no-memory 2: threads start with same seed and are reseeded thalfway with different seed, check final values are different."},
  {NULL,NULL,NULL}
 };

Ret_t testme( const tt& t) 
{
  G4Timer timer;
  timer.Start();
    std::cout<<std::left<<"\t  Sub-Test: ";
    std::cout<<std::left<<std::setfill(' ')<<std::setw(20)<<t.name;
  Ret_t result = t.aTest();
  if ( result == FAIL ) 
    {
      std::cout<<std::right<<" FAILED"<<std::flush;
      std::cerr<<"Sub-Test: "<<t.name<<":"<<t.mess<<" FAILED"<<std::endl;
    }
  else 
    {
        std::cout<<std::right<<" SUCCESS";
    }
  timer.Stop();
  std::cout<<"  (t="<<timer.GetUserElapsed()<<" s)"<<std::endl;
  return result;
}


int main(int,char**) 
{
    bool result = SUCCESS;
    std::cout<<"Testing Geant4 RNG"<<std::endl;
    std::cout<<"Running with "<<NTHREADS<<" threads, number histories: "<<TRIALS<<std::endl;
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
    for ( int t = 0 ; t < 5 ; ++t ) //Loop on engine types
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
