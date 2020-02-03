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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#include "G4MagneticField.hh"

#include <fstream>
#include <vector>

class TabulatedField3D

#ifndef STANDALONE
 : public G4MagneticField
#endif

{
    
public:
  
  TabulatedField3D(G4float gr1, G4float gr2, G4float gr3, G4float gr4, G4int quadModel);
  
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;

private:
    
  std::vector< std::vector< std::vector< double > > > fXField;
  
  std::vector< std::vector< std::vector< double > > > fYField;
  
  std::vector< std::vector< std::vector< double > > > fZField;
  
  G4int fNx,fNy,fNz; 
  
  G4double fMinix, fMaxix, fMiniy, fMaxiy, fMiniz, fMaxiz;
  
  G4double fDx, fDy, fDz;
  
  G4float fGradient1, fGradient2, fGradient3, fGradient4;
  
  G4int fModel;

};

