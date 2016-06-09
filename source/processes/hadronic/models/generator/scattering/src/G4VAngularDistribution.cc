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
// $Id: G4VAngularDistribution.cc,v 1.3 2003/02/05 15:43:42 gcosmo Exp $ //
//
// 
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4VAngularDistribution
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 15 May 1999
//
//      Modifications: 
//
// Id: G4VCrossSectionSource.cc,v 1.16 2000/05/11 19:07:29 pia Exp $ //   
//  
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4VAngularDistribution.hh"
#include "Randomize.hh"

G4double G4VAngularDistribution::Phi() const
{
  G4double phi = 2. * pi * G4UniformRand();
  return phi;
}
