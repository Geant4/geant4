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
// $Id: testReferenceCountedHandle.cc,v 1.5 2001-07-11 10:01:01 gunter Exp $
// 
//
// The program testing features and behaviour of the reference counted handle, smart-pointer, class.
// ----------------------------------------------------------------------
#include "globals.hh"
#include "G4ios.hh"
#include "G4ReferenceCountedHandle.hh"
#include <string>
#include <vector>

class TesterBase
{
public:  
  TesterBase() {}
  virtual ~TesterBase() {}
  bool Ok() const { return true; }
  virtual void Report() = 0;
};

static G4std::string sDefault = "Default";

class TesterString : public TesterBase
{
public:
  TesterString( G4std::string& str = sDefault ) : fData(str) {
  }
  
  ~TesterString(){
     G4cout << "Tester " << fData << " being destroyed..." << G4endl;    
  }
  
  const G4std::string& Data() const { return fData; }
  void  SetData( const std::string& data ) { fData = data; }
  void  SetData( const char*        data ) { fData = data; }
  virtual void Report() {
    G4cout << ">>" << Data() << "<< ";
  }
    
private:
  G4std::string fData;
};


class TesterInt : public TesterBase
{
public:
  TesterInt( int i = 0 ) : fData( i ) {
  }
  
  ~TesterInt(){
     G4cout << "Tester " << fData << " being destroyed..." << G4endl;    
  }
  
  const int& Data() const               { return fData; }
  void       SetData( const int data ) { fData = data; }
    
  virtual void Report() {
    G4cout << ">>" << Data() << "<< ";
  }

private:
  int fData;
};

typedef G4ReferenceCountedHandle<TesterBase>   Counted;
typedef G4ReferenceCountedHandle<TesterString> CountedString;
typedef G4ReferenceCountedHandle<TesterInt>    CountedInt;

//G4Allocator<G4ReferenceCountedHandle<TesterBase>::CountedObject>   G4ReferenceCountedHandle<TesterBase>::aRCHCountedObjectAllocator;
//G4Allocator<G4ReferenceCountedHandle<TesterString>::CountedObject> G4ReferenceCountedHandle<TesterString>::aRCHCountedObjectAllocator;
//G4Allocator<G4ReferenceCountedHandle<TesterInt>::CountedObject>    G4ReferenceCountedHandle<TesterInt>::aRCHCountedObjectAllocator;
//G4Allocator<G4ReferenceCountedHandle<TesterBase> >                 G4ReferenceCountedHandle<TesterBase>::aRCHAllocator;
//G4Allocator<G4ReferenceCountedHandle<TesterString> >               G4ReferenceCountedHandle<TesterString>::aRCHAllocator;
//G4Allocator<G4ReferenceCountedHandle<TesterInt> >                  G4ReferenceCountedHandle<TesterInt>::aRCHAllocator;

class HandleWatcher {
public:
  typedef std::vector<G4ReferenceCountedHandle<TesterBase> >       Handles;
  typedef Handles::iterator          HandlesIt;
  typedef Handles::const_iterator    HandlesCIt;

public:
  HandleWatcher() {
  }
  ~HandleWatcher() {
    fHandles.clear();
  }
  void AddHandle( const Counted& handle ) {
    SetCurrent( handle );
    fHandles.push_back( handle );
  }
  const Counted& operator[]( unsigned int index ) {
    return fHandles[index];
  }
  void SetCurrent( const Counted& handle ) {
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

int main( int argc, char* argv[] )
{
  G4cout << "Testing reference counting..." << G4endl;
  //HandleWatcher watcher;
  
  Counted t1 = new TesterString();
  {
    assert( t1 );
    assert( 1 == t1.Count() );

    G4cout << "Initial                 tester 1 count: " << t1.Count() << G4endl;
    
    //watcher.AddHandle( t1 );

    Counted t2;
    assert( !t2 );
    
    t2 = t1;

    assert( t1 );
    assert( t2 );
    //assert( 2 == t1.Count() );
    //assert( 2 == t2.Count() );

    G4cout << "After assignment        tester 1 count: " << t1.Count() << G4endl;
    G4cout << "Initial                 tester 2 count: " << t2.Count() << G4endl;

    //watcher.AddHandle( t2 );
    
    {
      Counted t3 = t2;

      assert( t3 );

      //assert( 3 == t1.Count() );
      //assert( 3 == t2.Count() );
      //assert( 3 == t3.Count() );

      G4cout << "After copy              tester 1 count: " << t1.Count() << G4endl;
      G4cout << "After copy              tester 2 count: " << t2.Count() << G4endl;
      G4cout << "Initial                 tester 3 count: " << t3.Count() << G4endl;
      
      //watcher.AddHandle( t3 );
    } // t3 Done

    assert( t1 );
    assert( t2 );
    //assert( 2 == t1.Count() );
    //assert( 2 == t2.Count() );

    G4cout << "After destruction of t3 tester 1 count: " << t1.Count() << G4endl;
    G4cout << "                        tester 2 count: " << t2.Count() << G4endl;
  } // t2 Done    
  
  assert( t1 );
  //assert( 1 == t1.Count() );
  
  G4cout << "After destruction of t2 tester 1 count: " << t1.Count() << G4endl;
  
  //--------------------------------------------------------------------------------
  G4cout << "Testing dereferencing..." << G4endl;
  assert( t1->Ok() );
  G4cout << "Dereferencing of t1   method Ok() is: " << (int)t1->Ok() << G4endl;
  
  const TesterBase* nativeObj = t1();
  assert( nativeObj != 0);
  G4cout << "Dereferenced object's method Ok() is: " << (int)nativeObj->Ok() << G4endl;
  
  //--------------------------------------------------------------------------------
  G4cout << "Testing copy and assignment by pointer to counted object..." << G4endl;
  G4cout << "Tester 1 is counting: >>"; t1->Report();
  G4cout << "<< with the count: " << t1.Count() << G4endl;
  Counted t4 = t1;
  G4cout << "Tester 1 is counting: >>"; t1->Report();
  G4cout << "<< with the count: " << t1.Count() << G4endl;
  G4cout << "Tester 4 is counting: >>"; t4->Report();
  G4cout << "<< with the count: " << t4.Count() << G4endl;

  //watcher.AddHandle( t4 ); 
  
  G4std::string str = "New one";
  t1 = new TesterString( str );
  assert( 1 == t1.Count() );    
  G4cout << "After assignment of native pointer to t1: " << G4endl;
  G4cout << "Tester 1 is counting: >>"; t1->Report();
  G4cout << "<< with the count: " << t1.Count() << G4endl;
  G4cout << "Tester 4 is counting: >>"; t4->Report();
  G4cout << "<< with the count: " << t4.Count() << G4endl;
  
  TesterString* dumbPtr = dynamic_cast<TesterString*>( t4() );
  dumbPtr->SetData( "Counted object's data updated" );
  G4cout << "Tester 4 is counting: >>"; t4->Report();
  G4cout << "<< with the count: " << t4.Count() << G4endl;
  t4 = new TesterInt( 100 );
  G4cout << "Tester 4 is counting: >>"; t4->Report();
  G4cout << "<< with the count: " << t4.Count() << G4endl;
  
  //--------------------------------------------------------------------------------
  
  //watcher.Report();
  
  return 0;
}
