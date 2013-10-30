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
// Unit test for G4AutoDelete
//
#include <iostream>
#include <cstdlib>
#include <map>
using std::map;
using std::cout;
using std::cerr;
using std::endl;

//#include "G4ThreadLocalSingleton.hh"
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

struct B {
  B() : a(-1) { MSG("B::B() Default constructor : "<<a); }
  B(int _a) : a(_a) { MSG("B::B(int) Constructor : "<<a); }
  ~B() { MSG("B::~B() Destructor : "<<a); }
  int a;
};

struct C {
  C() : a(-1) { MSG("C::C() Default constructor : "<<a); }
  C(int _a) : a(_a) { MSG("C::C(int) Constructor : "<<a); }
  ~C() { MSG("C::~C() Destructor : "<<a); }
  int a;
};

G4Mutex mapMutex = G4MUTEX_INITIALIZER;
map<pid_t,A*> amap;


G4ThreadFunReturnType myfunc( G4ThreadFunArgType val) {
  MSG("Starting thread");
  //Something here...
  A* a = new A( G4Threading::G4GetPidId() );
  G4AutoDelete::Register(a);
  G4AutoDelete::Register( new B(a->a) );
  G4AutoDelete::Register( new C(a->a) );
  G4AutoLock l(&mapMutex);
  amap[G4Threading::G4GetPidId()]=a;
  return NULL;
}


int main(int,char**) {
  MSG("Start of main function.");// Construction of variables at global scope has preceeded");
  //
 
  int nthreads = 2; //num threads
  G4Thread* tid = new G4Thread[nthreads];

  for ( int idx = 0 ; idx < nthreads ; ++idx ) {
    G4THREADCREATE( &(tid[idx]) , myfunc, static_cast<void*>(NULL) );
  }

  for ( int idx = 0 ; idx < nthreads ; ++idx ) {
    G4THREADJOIN( (tid[idx]) );
  }

  map<pid_t,A*>::iterator it = amap.begin();
  map<pid_t,A*>::iterator it2 = amap.begin();
  for ( ; it!=amap.end() ; ++it ) {
    cout<<it->first<<" "<<it->second<<" "<<it->second->a<<endl;
    CHECKVAL( it->first == it->second->a , "Wrong content of AutoDelete object");
  }

  MSG("End of main function, Exit, destruction of statics should follow");
  return 0;
}
