// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsTable.hh,v 1.6 2001-02-24 05:37:52 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//
//      Modified 01 March 1996, K. Amako
// ------------------------------------------------------------
//      STL migration Feb 2001, H.Kurashige

#ifndef G4PhysicsTable_h
#define G4PhysicsTable_h 1

#include "g4std/vector"
#include "globals.hh"
class     G4PhysicsVector;


class G4PhysicsTable: public G4std::vector<G4PhysicsVector*> 
{
 public: // with description

  G4PhysicsTable();
  // Deafult constructor.

  G4PhysicsTable(size_t capacity);
  // Constructor with capacity 


  virtual ~G4PhysicsTable();
  // Destructor

 private:
 // Private copy constructor and assignment operator.
  G4PhysicsTable(const G4PhysicsTable&);
  // Copy constructor  

  G4PhysicsTable& operator=(const G4PhysicsTable&);
  // assignment operator

 public:
  G4PhysicsVector*& operator()(size_t);
  G4PhysicsVector* const& operator()(size_t) const;

  void clearAndDestroy();
  // Remove all items and delete them at the same time

  void   insert (G4PhysicsVector*);
  void   insertAt (size_t, G4PhysicsVector*); 

  size_t entries() const;
  size_t length() const;
  G4bool isEmpty() const;
 
};

typedef G4PhysicsTable::iterator G4PhysicsTableIterator;

#include "G4PhysicsVector.hh"
#include "G4PhysicsTable.icc"
#endif
