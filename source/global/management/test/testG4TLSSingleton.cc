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
// $Id$
//
// ---------------------------------------------------------------
// Unit test for G4ThreadLocalSingleton

#include <iostream>
#include <cstdlib>
using std::cout;
using std::cerr;
using std::endl;

#include "G4ThreadLocalSingleton.hh"
#include "G4AutoDelete.hh"


#define CHECKVAL( cond , errmsg ) { \
    if( (cond) != true ) { \
      cerr<< errmsg << endl; \
      exit(1); \
    }\
}

G4Mutex aMutex = G4MUTEX_INITIALIZER;

#define MSG( msg ) { \
  G4AutoLock l(&aMutex); \
  cout<<G4Threading::G4GetPidId()<<" "<<msg<<endl;	\
  }


struct A {
  A() : a(-1) { MSG("A::A() Default constructor : "<<a); }
  A(int _a) : a(_a) { MSG("A::A(int) Constructor : "<<a); }
  ~A() { MSG("A::~A() Destructor : "<<a); }
  int a;
};

class G4SingletonExample {
  //ADDED
  friend class  G4ThreadLocalSingleton<G4SingletonExample>; 
private:
  //REMOVED
  //static G4SingletonExample* theInstance;
  G4SingletonExample() { MSG("This is G4SingletonExample constructor: "<<this); }
  static G4Mutex ctrm;
  static int ctr;
public:
  ~G4SingletonExample() { MSG("This is G4SingletonExample destructor: "<<this); }
  static G4SingletonExample* GetInstance() {
    //REPLACED
    static G4ThreadLocalSingleton<G4SingletonExample> inst;
    //static G4SingletonExample inst;
    // ADDED
    return inst.Instance();
    //REMOVED
    //if ( theInstance == 0 ) theInstance = &inst;
    //return theInstance;
  }
};
//REMOVED
//G4SingletonExample* G4SingletonExample::theInstance = 0;

G4ThreadFunReturnType myfunc( G4ThreadFunArgType val) {
  MSG("Starting thread");
  G4SingletonExample* inst = G4SingletonExample::GetInstance();
  MSG("In thread: "<<inst);
  CHECKVAL( inst == G4SingletonExample::GetInstance() , "Second call to G4SingletonExample::GetInstance() returns different address");
  //Something here...
  A* a = new A( G4Threading::G4GetPidId() );
  G4AutoDelete::Register(a);
  return NULL;
}


void foo(G4SingletonExample* test) {
  CHECKVAL( test == G4SingletonExample::GetInstance() , "Second call to G4SingletonExample::GetInstance() returns different address in function");
}

int main(int,char**) {
  MSG("Start of main function.");// Construction of variables at global scope has preceeded");
  //
  MSG("Testing G4ThreadLocalSingleton class, for class of type A");  
  static G4ThreadLocalSingleton<A> aSing;
  A* theS = aSing.Instance();
  CHECKVAL( theS == aSing.Instance() , "Second call to singleton returns different address");
  //
  MSG("Testing a singleton example");
  G4SingletonExample* inst = G4SingletonExample::GetInstance();
  CHECKVAL( inst == G4SingletonExample::GetInstance() , "Second call to G4SingletonExample::GetInstance() returns different address");
  foo(inst);

  int nthreads = 2; //num threads
  G4Thread* tid = new G4Thread[nthreads];

  for ( int idx = 0 ; idx < nthreads ; ++idx ) {
    G4THREADCREATE( &(tid[idx]) , myfunc, static_cast<void*>(inst) );
  }

  for ( int idx = 0 ; idx < nthreads ; ++idx ) {
    G4THREADJOIN( (tid[idx]) );
  }
  MSG("End of main function, Exit, destruction of statics should follow");
  return 0;
}
