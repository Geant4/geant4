// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testReferenceCountedHandle.cc,v 1.2 2001-04-18 19:40:50 radoone Exp $
// 
//
// The program testing features and behaviour of the reference counted handle, smart-pointer, class.
// ----------------------------------------------------------------------
#include "globals.hh"
#include "G4ios.hh"
#include "G4ReferenceCountedHandle.hh"
#include <string>

static G4std::string sDefault = "Default";

class Tester
{
public:
  Tester( G4std::string& str = sDefault ) : fData(str) {}
  ~Tester(){}
  bool Ok() { return true; }
  G4std::string& Data() { return fData; }
private:
  G4std::string fData;
};

typedef G4ReferenceCountedHandle<Tester> CountedTester;

int main( int argc, char* argv[] )
{
  G4cout << "Testing reference counting..." << G4endl;
  
  CountedTester t1 = new Tester();
  {
    assert( t1 );
    assert( 1 == t1.Count() );

    G4cout << "Initial                 tester 1 count: " << t1.Count() << G4endl;

    CountedTester t2;
    assert( !t2 );
    
    t2 = t1;

    assert( t1 );
    assert( t2 );
    assert( 2 == t1.Count() );
    assert( 2 == t2.Count() );

    G4cout << "After assignment        tester 1 count: " << t1.Count() << G4endl;
    G4cout << "Initial                 tester 2 count: " << t2.Count() << G4endl;

    {
      CountedTester t3 = t2;

      assert( t3 );

      assert( 3 == t1.Count() );
      assert( 3 == t2.Count() );
      assert( 3 == t3.Count() );

      G4cout << "After copy              tester 1 count: " << t1.Count() << G4endl;
      G4cout << "After copy              tester 2 count: " << t2.Count() << G4endl;
      G4cout << "Initial                 tester 3 count: " << t3.Count() << G4endl;
    } // t3 Done

    assert( t1 );
    assert( t2 );
    assert( 2 == t1.Count() );
    assert( 2 == t2.Count() );

    G4cout << "After destruction of t3 tester 1 count: " << t1.Count() << G4endl;
    G4cout << "                        tester 2 count: " << t2.Count() << G4endl;
  } // t2 Done    
  
  assert( t1 );
  assert( 1 == t1.Count() );
  
  G4cout << "After destruction of t2 tester 1 count: " << t1.Count() << G4endl;
  
  //--------------------------------------------------------------------------------
  G4cout << "Testing dereferencing..." << G4endl;
  assert( t1->Ok() );
  G4cout << "Dereferencing of t1   method Ok() is: " << (int)t1->Ok() << G4endl;
  
  Tester* nativeObj = t1();
  assert( nativeObj != 0);
  G4cout << "Dereferenced object's method Ok() is: " << (int)nativeObj->Ok() << G4endl;
  
  //--------------------------------------------------------------------------------
  G4cout << "Testing copy and assignment by pointer to counted object..." << G4endl;
  G4cout << "Tester 1 is counting: >>" << t1->Data() << "<< with the count: " << t1.Count() << G4endl;
  G4std::string str;
  t1 = new Tester( str = "New one" );
  assert( 1 == t1.Count() );    
  G4cout << "Tester 1 is counting: >>" << t1->Data() << "<< with the count: " << t1.Count() << G4endl;
  
  //--------------------------------------------------------------------------------
  
  return 0;
}
