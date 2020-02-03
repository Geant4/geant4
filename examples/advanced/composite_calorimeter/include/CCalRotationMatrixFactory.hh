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
///////////////////////////////////////////////////////////////////////////////
// File: CCalRotationMatrixFactory.hh
// Description: CCalRotationFactory is a singleton class to get from a file
//              the information to build CCal Rotation matrices.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalRotationMatrixFactory_h
#define CCalRotationMatrixFactory_h 1

#include <map>
#include "G4RotationMatrix.hh"

typedef G4RotationMatrix* G4RotationMatrixPtr;
typedef std::map<G4String, G4RotationMatrixPtr, std::less<G4String> > G4RotationMatrixTable;
typedef std::map<G4String, G4RotationMatrixPtr, std::less<G4String> >::iterator G4RotationMatrixTableIterator;


//typedef RWTPtrOrderedVector<G4RotationMatrix> G4RotationMatrixTable;
//  Is an instantiation of the template class RWTPtrOrderedVector.

class CCalRotationMatrixFactory {
public:
  ~CCalRotationMatrixFactory();

  static CCalRotationMatrixFactory* getInstance();
  static CCalRotationMatrixFactory* getInstance(const G4String& rotfile);
  static void setFileName(const G4String& rotfile);

  G4RotationMatrix* findMatrix(const G4String&);
  G4RotationMatrix* AddMatrix(const G4String& name, 
                              G4double th1, G4double phi1,  //Axis angles
                              G4double th2, G4double phi2,  //in rads
                              G4double th3, G4double phi3); //

private:
  CCalRotationMatrixFactory();

private:
  static CCalRotationMatrixFactory* instance;
  static G4String file;

  G4RotationMatrixTable theMatrices; //Where the matrices are stored.
};

// 29-Jan-2004 A.R. : commented to avoid clashes with CLHEP.
//                    Streaming operators for rotation matrices are
//                    already defined in CLHEP::HepRotation.
// std::ostream& operator<<(std::ostream&, const G4RotationMatrix &);
#endif
