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
// $Id: G4DataVector.cc,v 1.8 2006/06/29 19:03:59 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//  G4DataVector.cc
//
//  History:
//    18 Sep. 2001, H.Kurashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4DataVector.hh"
#include <iomanip>

G4DataVector::G4DataVector()
  : std::vector<G4double>()
{
}

G4DataVector::G4DataVector(size_t capacity)
  : std::vector<G4double>(capacity, 0.0)
{
}

G4DataVector::G4DataVector(size_t capacity, G4double value)
  : std::vector<G4double>(capacity, value)
{
}

G4DataVector::~G4DataVector()
{
}

G4bool G4DataVector::Store(std::ofstream& fOut, G4bool ascii)
{
  // Ascii mode
  if (ascii)
  {
    fOut << *this;
    return true;
  } 

  // Binary Mode
  size_t sizeV = size(); 
  fOut.write((char*)(&sizeV), sizeof sizeV);

  G4double* value = new G4double[sizeV];
  size_t i=0;
  for (const_iterator itr=begin(); itr!=end(); itr++, i++)
  {
    value[i]  =  *itr;
  }
  fOut.write((char*)(value), sizeV*(sizeof (G4double)) );
  delete [] value;
  
  return true;
}

G4bool G4DataVector::Retrieve(std::ifstream& fIn, G4bool ascii)
{
  clear();
  size_t sizeV=0;
  
  // retrieve in ascii mode
  if (ascii)
  {
    // contents
    fIn >> sizeV;
    if (fIn.fail())  { return false; }
    
    reserve(sizeV);
    for(size_t i = 0; i < sizeV ; i++)
    {
      G4double vData=0.0;
      fIn >> vData;
      if (fIn.fail())  { return false; }
      push_back(vData);
    }
    return true ;
  }
  
  // retrieve in binary mode
  fIn.read((char*)(&sizeV), sizeof sizeV); 
 
  G4double* value = new G4double[sizeV];
  fIn.read((char*)(value),  sizeV*(sizeof(G4double)) );
  if (G4int(fIn.gcount()) != G4int(sizeV*(sizeof(G4double))) )
  {
    delete [] value;
    return false;
  }

  reserve(sizeV);
  for(size_t i = 0; i < sizeV; i++) 
  {
    push_back(value[i]);
  }
  delete [] value;
  return true;
}
    
std::ostream& operator<<(std::ostream& out, const G4DataVector& pv)
{
  out << pv.size() << G4endl; 
  for(size_t i = 0; i < pv.size(); i++)
  {
    out << std::setprecision(12) << pv[i] << G4endl;
  }

  return out;
}

