// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DataVector.hh,v 1.7 2001-03-06 15:56:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
// ------------------------------------------------------------
//
// Class Description:
//
//   Utility class providing similar behaviour of vector<G4double>.
//   It includes additional methods for compatibility with Rogue Wave
//   collection.
//

#ifndef G4DataVector_h
#define G4DataVector_h 1

#include "globals.hh"
#include "g4std/vector"

class G4DataVector : public G4std::vector<G4double> 
{

 public: // with description

  G4DataVector();
   // Default constructor.

  G4DataVector(size_t capacity);
   // Constructor given a 'capacity' defining the initial number of elements.

  virtual ~G4DataVector(){;}
   // Empty destructor

  inline void insertAt(size_t, const G4double&);
    // Insert an element at given position

  inline size_t index(const G4double&);
    // Returns back index of the element same as given value

  inline G4bool contains(const G4double&) const;
    // Returns 'true' if it contains the element same as given value 

  inline G4bool remove(const G4double&);
    // Removes the first element same as given value  

  inline size_t removeAll(const G4double&);
    // Remove all elements same as given value  
};

inline
G4DataVector::G4DataVector()
  : G4std::vector<G4double>()
{
}

inline
G4DataVector::G4DataVector(size_t capacity)
  : G4std::vector<G4double>(capacity)
{
}

inline
void G4DataVector::insertAt(size_t pos, const G4double& a) 
{ 
  iterator i = begin();
  for (size_t ptn=0; (ptn<pos)&&(i!=end()); i++,ptn++);
  if (i==end())
    push_back(a);
  else
    insert(i,a);
}

inline
size_t G4DataVector::index(const G4double& a) 
{ 
  size_t ptn = 0;
  for (iterator i=begin(); i!=end(); i++,ptn++)
    if (*i==a) return ptn;

  return (ptn=~(size_t)0);
}

inline
G4bool G4DataVector::contains(const G4double& a) const
{
  for (const_iterator i=begin(); i!=end(); i++)
    if (*i==a) return true;

  return false;
}

inline 
G4bool G4DataVector::remove(const G4double& a)
{
  G4bool found =false;

  for (iterator i=begin(); i!=end(); i++){
    if (*i==a) {
      erase(i);
      found = true;
      break;
    }
  }
  return found;
}

inline 
size_t G4DataVector::removeAll(const G4double& a)
{
  size_t ptn=0;

  for (iterator i=begin(); i!=end(); i++){
    if (*i==a) {
      erase(i);
      ptn++;
      i--;
    } 
  }
  return ptn;
}

#endif
