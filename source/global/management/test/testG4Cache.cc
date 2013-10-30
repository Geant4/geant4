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
// Unit test for G4Cache
#include <iostream>
#include <cstdlib>
using std::cout;
using std::cerr;
using std::endl;

#include "G4Cache.hh"

struct A {
  A(int _a) : a(_a) {}
  int a;
};

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

G4ThreadFunReturnType myfunc( G4ThreadFunArgType val) {
  MSG("Starting thread");
  G4MapCache<int,double>& asc = *((G4MapCache<int,double>*)val);
  asc[1]=G4Threading::G4GetPidId()%100000;
  asc[2]=2.2;
  asc[30]=3.3;
  CHECKVAL( asc[1]==G4Threading::G4GetPidId()%100000,"Wrong first element");
  CHECKVAL( asc[2]==2.2,"Wrong second element");
  CHECKVAL( asc[30]==3.3,"Wrong third element");
  MSG("Shared Cache object content from thread:");
  for ( G4MapCache<int,double>::const_iterator it = asc.Begin() ; it!=asc.End() ; ++it ) 
    MSG(it->first<<":"<<it->second);
  MSG("Thread is over");
  return NULL;
}

int main(int,char**) {
  MSG("Testing G4caches");
  int nthreads = 0; //num threads
  G4Thread* tid = new G4Thread[nthreads];

  G4MapCache<int,double> aSharedCache;
  aSharedCache.Insert(1,0.1);
  aSharedCache.Insert(2,0.2);
  for ( int idx = 0 ; idx < nthreads ; ++idx ) {
    G4THREADCREATE( &(tid[idx]) , myfunc, &(aSharedCache) );
  }

  for ( int idx = 0 ; idx < nthreads ; ++idx ) {
    G4THREADJOIN( (tid[idx]) );
  }

  cout<<"Shared Cache object content from main:"<<endl;
  for ( G4MapCache<int,double>::const_iterator it = aSharedCache.Begin() ; it!=aSharedCache.End() ; ++it ) cout<<it->first<<":"<<it->second<<endl;

  cout<<"Test single value cache"<<endl;
  G4Cache<double> vc(3.431);
  CHECKVAL( vc.Get() == 3.431 , "Wrong content of cache");
  double theV = vc.Get();
  theV += 1.5;
  vc.Put( theV );
  CHECKVAL( vc.Get() == 4.931 , "Wrong manipulation of cache");
  cout<<"Final conten of cache"<<endl;
  cout<<vc.Get()<<endl;
  
  cout<<"Test vector cache"<<endl;
  G4VectorCache<double> aV;
   aV.Push_back( 1.01 );
  CHECKVAL( aV.Size() == 1 , "Wrong size, expected 1" );
  CHECKVAL( aV[0]==1.01, "Wrong content for first element");
  CHECKVAL( aV.Pop_back()==1.01,"Wrong content of head");
  CHECKVAL( aV.Size() == 0 , "Wrong size, expected 0" );
  aV.Push_back(3);
  aV[0]=2;
  CHECKVAL( aV[0]==2,"Wrong replacement of value");
  cout<<"Array constructor"<<endl;
  double array[3]={1.1,2.2,3.3};
  G4VectorCache<double> aV2(3,array);
  CHECKVAL( (aV2[0]==1.1)&&(aV2[1]==2.2)&&(aV2[2]==3.3),"Wrong content of array constructor");
  cout<<"Content of array cache:"<<endl;
  cout<<aV2[0]<<" "<<aV2[1]<<" "<<aV2[2]<<endl;
  //Check "stack" funcionality
  aV2.Push_back(4.4);
  CHECKVAL( aV2.Pop_back()==4.4,"Wrong last element");
  CHECKVAL( aV2.Pop_back()==3.3,"Wrong last element");
  CHECKVAL( aV2.Pop_back()==2.2,"Wrong last element");
  CHECKVAL( aV2.Pop_back()==1.1,"Wrong last element"); 

  cout<<"Test map cache"<<endl;
  G4MapCache<int,double> aM1,aM2;
  aM1[1]=10.1;
  aM1[10]=100.1;
  aM1[5]=50.1;
  CHECKVAL( aM1.Size()==3,"Wrong map size");
  std::pair<G4MapCache<int,double>::iterator,bool> e = aM1.Insert(20,200.1);
  CHECKVAL( (e.second==true)&&(e.first->first==20)&&(e.first->second==200.1),"Wrong insert");
  CHECKVAL( (aM1.Begin()->first==1)&&(aM1.Begin()->second==10.1),"Wrong head");
  CHECKVAL( aM1.Find(5)!=aM1.End() , "Wrong find");
  CHECKVAL( aM1.Find(111)==aM1.End() , "Wrong find");
  CHECKVAL( aM1.Has(5)==true,"Wrong Has");
  CHECKVAL( aM1.Get(10)==100.1,"Wrong Get");
  CHECKVAL( (aM1.Erase(10)==1)&&(aM1.Has(10)==false),"Wrong Erase");
  CHECKVAL( aM1[20]==200.1,"Wrong []");
  aM1[20]=199.9;
  CHECKVAL( aM1[20]==199.9,"Wrong []");
  cout<<"Content of map cache:"<<endl;
  for ( G4MapCache<int,double>::const_iterator it = aM1.Begin() ; it!=aM1.End() ; ++it )
    cout<<it->first<<":"<<it->second<<endl;

  //Try with class
  cout<<"Test map-cache with pointers to objects"<<endl;
  G4MapCache<int,A*> mm;
  mm[0]=new A(10);
  mm[1]=new A(11);
  mm[2]=new A(12);
  CHECKVAL( mm[0]->a==10 , "Wrong first element");
  CHECKVAL( mm[1]->a==11 , "Wrong second element");
  CHECKVAL( mm[2]->a==12 , "Wrong third element");
  cout<<"Content of map cache for objects"<<endl;
  for ( G4MapCache<int,A*>::const_iterator it = mm.Begin() ; it!=mm.End() ; ++it )
    cout<<it->first<<":"<<it->second<<"="<<it->second->a<<endl;
    
  //Try with class pointer!
  cout<<"Test simple cache with pointer to object"<<endl;
  G4Cache<A*> ac;
  A* mya = 0;
  ac.Put( mya=new A(2) );
  mya->a = 1;
  CHECKVAL( ac.Get()->a == 1 , "Wrong Cache content");
  
  cout<<"SUCCESS"<<endl;
  return 0;
}
