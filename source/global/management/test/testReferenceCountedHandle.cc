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
// $Id: testReferenceCountedHandle.cc,v 1.8 2001-11-07 14:13:58 radoone Exp $
// 
//
// The program testing features and behaviour of the
// reference counted handle, smart-pointer, class.
// ----------------------------------------------------------------------

#define sDefault  "Default"
#define csDefault "cDefault"

#include "G4ios.hh"
#include "G4ReferenceCountedHandle.hh"

#include <assert.h>

#include <iostream>
#include <string>
#include <vector>

class TesterBase
{
public:  
  TesterBase() {
  }
  
  virtual ~TesterBase() {
    G4cout << "TesterBase being destroyed..." << G4endl;
  }
  
  inline bool Ok() const {
    return true;
  }
  
  virtual void Report() = 0;
};

class TesterString : public TesterBase
{
public:
  
  TesterString( const char* str = csDefault ) {
    fData = str;
  }
  
  TesterString( G4std::string& str ) : fData(str) {}
  
  ~TesterString(){
    G4cout << "Tester " << fData << " being destroyed..." << G4endl;    
  }
  
  inline const G4std::string& Data() const {
    return fData;
  }
  
  inline void  SetData( const G4std::string& data ) {
    fData = data;
  }
  
  inline void  SetData( const char*        data ) {
    fData = data;
  }
  
  virtual void Report() {
    G4cout << " >>" << Data() << "<< ";
  }
  
private:
  G4std::string fData;
};


class TesterInt : public TesterBase
{
public:
  TesterInt( int i = 0 ) : fData( i ) {}
  
  ~TesterInt(){
    G4cout << "Tester " << fData << " being destroyed..." << G4endl;    
  }
  
  inline const int& Data() const {
    return fData;
  }
  
  inline void SetData( const int data ) {
    fData = data;
  }
  
  virtual void Report() {
    G4cout << " >>" << Data() << "<< ";
  }
  
private:
  int fData;
};

typedef G4ReferenceCountedHandle<TesterBase>   Counted;
typedef G4ReferenceCountedHandle<TesterString> CountedString;
typedef G4ReferenceCountedHandle<TesterInt>    CountedInt;

/*
// Beta prototype of the class monitoring run-time activity of the RC handles
class HandleWatcher {
public:
  typedef G4std::vector<G4ReferenceCountedHandle<TesterBase> >       Handles;
  typedef Handles::iterator          HandlesIt;
  typedef Handles::const_iterator    HandlesCIt;
  
public:
  HandleWatcher() {}
  
  ~HandleWatcher() {
    fHandles.clear();
  }
  
  inline void AddHandle( const Counted& handle ) {
    SetCurrent( handle );
    fHandles.push_back( handle );
  }
  
  inline const Counted& operator[]( unsigned int index ) {
    return fHandles[index];
  }
  
  inline void SetCurrent( const Counted& handle ) {
    fCurrent = handle;
  }
  
  void Report() {
    G4cout << "HandleWatcher Report:" << G4endl;
    G4cout << "---------------------" << G4endl;
    HandlesCIt it;
    for( it = fHandles.begin(); it != fHandles.end(); it++ ) {
      G4cout << "Handle "; (*it)->Report(); G4cout << " count " << (*it).Count() << G4endl;
    }
    G4cout << "\nCurrent handle "; fCurrent->Report();
    G4cout << " count " << fCurrent.Count() << G4endl;
  }
  
private:
  Handles fHandles;
  Counted fCurrent;
};
*/

inline void SimulateNavigator( Counted& rcObj ) {
  rcObj = new TesterString( "Simulated" );
}

int main( int argc, char* argv[] )
{
  Counted t0;
  G4cout << "Testing reference counting..." << G4endl;
  //HandleWatcher watcher;
  
  Counted t1 = new TesterString();
  t0 = t1;
  {
		  assert( t1 );
      assert( 2 == t1.Count() );
      
      G4cout << G4endl << "Initial tester 1 count:................... " << t1.Count() << G4endl;
      
      //watcher.AddHandle( t1 );
      
      Counted t2;
      assert( !t2 );
      
      G4cout << "Initial tester 2 count:................... " << t2.Count() << G4endl;
      
      t2 = t1;
      
      assert( t1 );
      assert( t2 );
      assert( 3 == t1.Count() );
      assert( 3 == t2.Count() );
      
      G4cout << "After assignment tester 1 count:.......... " << t1.Count() << G4endl;
      G4cout << "After assignment tester 2 count:.......... " << t2.Count() << G4endl;
      
      //memory_watcher.AllocInfo();
      //watcher.AddHandle( t2 );
      
      {
        Counted t3 = t2;
        
        assert( t3 );
        
        assert( 4 == t1.Count() );
        assert( 4 == t2.Count() );
        assert( 4 == t3.Count() );
        
        G4cout << "After copy tester 1 count:................ " << t1.Count() << G4endl;
        G4cout << "After copy tester 2 count:................ " << t2.Count() << G4endl;
        G4cout << "Initial tester 3 count:................... " << t3.Count() << G4endl;
        
        //watcher.AddHandle( t3 );
      } // t3 Done
      
      //memory_watcher.AllocInfo();
      
      assert( t1 );
      assert( t2 );
      assert( 3 == t1.Count() );
      assert( 3 == t2.Count() );
      
      G4cout << "After destruction of t3 tester 1 count:... " << t1.Count() << G4endl;
      G4cout << "                        tester 2 count:... " << t2.Count() << G4endl;
  } // t2 Done    
  
  assert( t1 );
  assert( 2 == t1.Count() );
  
  G4cout << "After destruction of t2 tester 1 count:... " << t1.Count() << G4endl;
  
  //--------------------------------------------------------------------------------
  G4cout << "\nTesting dereferencing..." << G4endl;
  assert( t1->Ok() );
  G4cout << "Dereferencing of t1 method Ok() is:....... " << (int)t1->Ok() << G4endl;
  
  const TesterBase* nativeObj = t1();
  assert( nativeObj != 0);
  G4cout << "Dereferenced object's method Ok() is:..... " << (int)nativeObj->Ok() << G4endl;
  
  //--------------------------------------------------------------------------------
  G4cout << "\nTesting copy and assignment\nby a pointer to counted object...\n" << G4endl;
  G4cout << "Tester 1 is counting:....................."; t1->Report();
  G4cout << "with count: " << t1.Count() << G4endl;
  Counted t4 = t1;
  G4cout << "\nTester 1 being copied to the new tester 4 ...\n" << G4endl;
  G4cout << "Tester 1 is counting:....................."; t1->Report();
  G4cout << "with count: " << t1.Count() << G4endl;
  G4cout << "Tester 4 is counting:....................."; t4->Report();
  G4cout << "with count: " << t4.Count() << G4endl;
  
  //watcher.AddHandle( t4 ); 
  
  G4std::string str = "New";
  t1 = new TesterString( str );
  assert( 1 == t1.Count() );    
  G4cout << "\nAfter assignment of a counted object pointer to t1:\n" << G4endl;
  G4cout << "Tester 0 is counting:....................."; t0->Report();
  G4cout << "with count: " << t0.Count() << G4endl;
  G4cout << "Tester 1 is counting:....................."; t1->Report();
  G4cout << "with count: " << t1.Count() << G4endl;
  G4cout << "Tester 4 is counting:....................."; t4->Report();
  G4cout << "with count: " << t4.Count() << G4endl;
  
  TesterString* dumbPtr = dynamic_cast<TesterString*>( t4() );
  dumbPtr->SetData( "Updated" );
  G4cout << "\nAfter the counted object's data update in t4:\n" << G4endl;
  G4cout << "Tester 0 is counting:....................."; t0->Report();
  G4cout << "with count: " << t0.Count() << G4endl;
  G4cout << "Tester 1 is counting:....................."; t1->Report();
  G4cout << "with count: " << t1.Count() << G4endl;
  G4cout << "Tester 4 is counting:....................."; t4->Report();
  G4cout << "with count: " << t4.Count() << G4endl;
  
  t4 = new TesterInt( 100 );
  G4cout << "\nAfter assignment of an int tester object to t1:\n" << G4endl;
  G4cout << "Tester 0 is counting:....................."; t0->Report();
  G4cout << "with count: " << t0.Count() << G4endl;
  G4cout << "Tester 1 is counting:....................."; t1->Report();
  G4cout << "with count: " << t1.Count() << G4endl;
  G4cout << "Tester 4 is counting:....................."; t4->Report();
  G4cout << "with count: " << t4.Count() << G4endl;
  
  G4cout << "\nCalling SimulateNavigator( t0 )...\n" << G4endl;
  SimulateNavigator( t0 );
  
  G4cout << "\nAfter SimulateNavigator( t0 ) call...\n" << G4endl;
  G4cout << "Tester 0 is counting:....................."; t0->Report();
  G4cout << "with count: " << t0.Count() << G4endl;
  G4cout << "Tester 1 is counting:....................."; t1->Report();
  G4cout << "with count: " << t1.Count() << G4endl;
  G4cout << "Tester 4 is counting:....................."; t4->Report();
  G4cout << "with count: " << t4.Count() << G4endl;
  
  //watcher.Report();
  
  G4cout << "\nEnd of program, expecting the automatic heap memory cleanup...\n" << G4endl;
  
  return 0;
}
