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
// G4GDMLReadDefine
//
// Class description:
//
// GDML class for positionings and transformations according to
// specifications in Geant4.

// Author: Zoltan Torzsok, November 2007
// --------------------------------------------------------------------
#ifndef G4GDMLREADDEFINE_HH
#define G4GDMLREADDEFINE_HH 1

#include <map>

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4GDMLRead.hh"

class G4GDMLMatrix
{
  public:

    G4GDMLMatrix();
    G4GDMLMatrix(std::size_t rows0, std::size_t cols0);
    G4GDMLMatrix(const G4GDMLMatrix& rhs);
    G4GDMLMatrix& operator=(const G4GDMLMatrix& rhs);
    ~G4GDMLMatrix();

    void Set(std::size_t r, std::size_t c, G4double a);
    G4double Get(std::size_t r, std::size_t c) const;
    std::size_t GetRows() const;
    std::size_t GetCols() const;

  private:

    G4double* m = nullptr;
    std::size_t rows = 0, cols = 0;
};

class G4GDMLReadDefine : public G4GDMLRead
{
  public:

    G4bool IsValidID(const G4String&) const;
    G4double GetConstant(const G4String&);
    G4double GetVariable(const G4String&);
    G4double GetQuantity(const G4String&);
    G4ThreeVector GetPosition(const G4String&);
    G4ThreeVector GetRotation(const G4String&);
    G4ThreeVector GetScale(const G4String&);
    G4GDMLMatrix GetMatrix(const G4String&);

    virtual void DefineRead(const xercesc::DOMElement* const);

    inline void SetReverseSearch(G4bool flag)
    { reverseSearch = flag; }

  protected:

    G4GDMLReadDefine();
    virtual ~G4GDMLReadDefine();

    G4RotationMatrix GetRotationMatrix(const G4ThreeVector&);
    void VectorRead(const xercesc::DOMElement* const, G4ThreeVector&);
    G4String RefRead(const xercesc::DOMElement* const);

    void ConstantRead(const xercesc::DOMElement* const);
    void MatrixRead(const xercesc::DOMElement* const);
    void PositionRead(const xercesc::DOMElement* const);
    void RotationRead(const xercesc::DOMElement* const);
    void ScaleRead(const xercesc::DOMElement* const);
    void VariableRead(const xercesc::DOMElement* const);
    void QuantityRead(const xercesc::DOMElement* const);
    void ExpressionRead(const xercesc::DOMElement* const);

  protected:
    G4bool reverseSearch = false;
    std::map<G4String, G4double> quantityMap;
    std::map<G4String, G4ThreeVector> positionMap;
    std::map<G4String, G4ThreeVector> rotationMap;
    std::map<G4String, G4ThreeVector> scaleMap;
    std::map<G4String, G4GDMLMatrix> matrixMap;
};

#endif
