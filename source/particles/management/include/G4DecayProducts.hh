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
// $Id: G4DecayProducts.hh,v 1.8 2001-07-11 10:01:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      12 Dec 1997 H.Kurashige
// ------------------------------------------------------------

#ifndef G4DecayProducts_h
#define G4DecayProducts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4Allocator.hh"

class G4DecayProducts
{
  public: // With Description

    // constructors
    G4DecayProducts();
    G4DecayProducts(const G4DynamicParticle &aParticle);

  public: 
    // copy constructor and assignment operator 
    //   Deep    copy: for G4DynamicParticle pointer
    G4DecayProducts(const G4DecayProducts &right);
    G4DecayProducts & operator=(const G4DecayProducts &right);

    //destructor
    ~G4DecayProducts();

    // new & delete operators for G4Allocator
    inline void *operator new(size_t);
    inline void operator delete(void *G4DecayProducts);

    // (un)equal operator
    G4int operator==(const G4DecayProducts &right) const;
    G4int operator!=(const G4DecayProducts &right) const;

  public: // With Description
   //  set-get methods for the parent particle   
   //    theParentPaticle is used to get information of parent particle 
   //    when decay products are filled 
   //    new G4DynamicParticle object is created in set methods  
    const G4DynamicParticle* GetParentParticle() const {return theParentParticle;};
    void SetParentParticle(const G4DynamicParticle &aParticle);

   //  boost all products
    void Boost(G4double totalEnergy, const G4ThreeVector &momentumDirection);
    void Boost(G4double betax, G4double betay, G4double betaz);
 
  //   push-pop  methods for decay products pointer
    G4DynamicParticle* PopProducts();
    G4int PushProducts(G4DynamicParticle *aParticle);

    G4DynamicParticle* operator[](G4int anIndex) const;

    G4int entries() const {return numberOfProducts;};

  // check energy/momentum of products 
    G4bool IsChecked() const; 
   
  // 
    void DumpInfo() const;

  protected:
    enum {MaxNumberOfProducts = 64};

  private: 
    G4int                         numberOfProducts;
    G4DynamicParticle*            theParentParticle;
    G4DynamicParticle*            theProductVector[MaxNumberOfProducts];

};

// ------------------------
// Inlined operators
// ------------------------
extern G4Allocator<G4DecayProducts> aDecayProductsAllocator;

inline void * G4DecayProducts::operator new(size_t)
{
  void * aDecayProducts;
  aDecayProducts = (void *) aDecayProductsAllocator.MallocSingle();
  return aDecayProducts;
}

inline void G4DecayProducts::operator delete(void * aDecayProducts)
{
  aDecayProductsAllocator.FreeSingle((G4DecayProducts *)  aDecayProducts);
}

inline 
 G4int G4DecayProducts::operator==(const G4DecayProducts &right) const
{
  return (this == (G4DecayProducts *) &right);
}

inline 
 G4int G4DecayProducts::operator!=(const G4DecayProducts &right) const
{
  return (this != (G4DecayProducts *) &right);
}

#endif

