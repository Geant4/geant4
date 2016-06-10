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
//
// $Id: G4DecayProducts.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      12 Dec 1997 H.Kurashige
//
//      Use std::vector 4  Apr. 2012  H.Kurashige
// ------------------------------------------------------------

#ifndef G4DecayProducts_h
#define G4DecayProducts_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include  <vector>

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

    typedef std::vector<G4DynamicParticle*>  G4DecayProductVector;
  protected:
    enum {MaxNumberOfProducts = 64};

  private: 
    G4int                               numberOfProducts;
    G4DynamicParticle*                  theParentParticle;
    G4DecayProductVector*               theProductVector;

};

// ------------------------
// Inlined operators
// ------------------------

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

