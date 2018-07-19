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
// $Id: G4PhysicsLinearVector.hh 98864 2016-08-15 11:53:26Z gcosmo $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class header file
//
//  G4PhysicsLinearVector.hh
//
//  Class description:
//
//    A physics vector which has values of energy-loss, cross-section, 
//    and other physics values of a particle in matter in a given 
//    range of the energy, momentum, etc. The scale of energy/momentum
//    bins is in linear.

//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    03 Mar. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Cache mechanism and hidden bin from the 
//                            user introduced
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added
//    11 Nov. 2000, H.Kurashige : Use STL vector for dataVector and binVector
//    16 Aug. 2011  H.Kurashige : Move dBin, baseBin to the base class
//    02 Oct. 2013  V.Ivanchenko : Remove FindBinLocation method
//
//--------------------------------------------------------------------

#ifndef G4PhysicsLinearVector_h
#define G4PhysicsLinearVector_h 1

#include "globals.hh"
#include "G4PhysicsVector.hh"

class G4PhysicsLinearVector : public G4PhysicsVector  
{
public:// with description

  G4PhysicsLinearVector();
       // the vector will be filled from external file using Retrieve method

  G4PhysicsLinearVector(G4double theEmin, G4double theEmax, size_t theNbin);
       // Energy vector will be computed and filled at construction, 
       // number of elements 'theNbin+1'. Use PutValue() to fill the data vector

  virtual ~G4PhysicsLinearVector();

  virtual G4bool Retrieve(std::ifstream& fIn, G4bool ascii) final;
       // To retrieve persistent data from a file stream.

  virtual void ScaleVector(G4double factorE, G4double factorV) final;
       // Scale all values of the vector and second derivatives
       // by factorV, energies - by vectorE. 

};


#endif
