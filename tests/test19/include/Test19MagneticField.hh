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
// $Id: Test19MagneticField.hh,v 1.1 2009-11-17 17:03:15 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Created by M. Kosov 11-Nov-2009 (following N02 Example)
//
//-------------------------------------------------------------------- 
//  A class for control of the Uniform Magnetic Field of the detector.
//-------------------------------------------------------------------- 


#ifndef Test19MagneticField_H
#define Test19MagneticField_H

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

class Test19MagneticField: public G4UniformMagField
{
  public:
  
   Test19MagneticField(G4ThreeVector);  //  The value of the field
   Test19MagneticField();               //  A zero field (dummy initialization)
  ~Test19MagneticField(){}  
      
   
   void SetMagFieldValue(G4double BValue) {SetMagFieldValue(G4ThreeVector(BValue,0,0));}
   void SetMagFieldValue(G4ThreeVector fieldVector);
      
   G4ThreeVector GetConstantFieldValue();

  protected:

   G4FieldManager* GetGlobalFieldManager()
           {return G4TransportationManager::GetTransportationManager()->GetFieldManager();}
};

#endif
