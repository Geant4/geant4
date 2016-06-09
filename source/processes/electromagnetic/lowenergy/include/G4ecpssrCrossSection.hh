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
// $Id: G4ecpssrCrossSection.hh,v 1.3 2008/06/03 07:29:15 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Author: Haifa Ben Abdelouahed
//         
//
// History:
// -----------
//  21 Apr 2008   H. Ben Abdelouahed   1st implementation
//  21 Apr 2008   MGP        Major revision according to a design iteration
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics, Cross section, p ionisation, K shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------


#ifndef G4ECPSSRCROSSSECTION_HH
#define G4ECPSSRCROSSSECTION_HH 1

#include "globals.hh"


class G4ecpssrCrossSection 

{
public:

  G4ecpssrCrossSection();

  virtual ~G4ecpssrCrossSection();
			     
  G4double CalculateCrossSection(G4int zTarget,G4int zIncident, G4double energyIncident);//according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)
				    
  G4double CalculateVelocity(G4int zTarget,G4int zIncident, G4double energyIncident); 
			      
  G4double  ExpIntFunction(G4int n,G4double x);//Exponential Integral Function

   

private:


  G4ecpssrCrossSection(const G4ecpssrCrossSection&);
  G4ecpssrCrossSection & operator = (const G4ecpssrCrossSection &right);

};

#endif
