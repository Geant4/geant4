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
// $Id: G4NystromRK4.hh,v 1.1 2009-11-05 11:29:25 japost Exp $ 
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NystromRK4
//
// Class description:
//
// Integrate the equations of the motion of a particle in a magnetic field
// using 4th Runge-Kutta-Nystrom method with errors estimation 
// (ATL-SOFT-PUB-2009-01)
// Current form can be used only for 'pure' magnetic field.
// 
// History:
// - Created: I.Gavrilenko   15.05.2009   (as G4AtlasRK4)
// - Adaptations:  J. Apostolakis  May-Nov 2009
// -------------------------------------------------------------------

#ifndef G4ATLASRK4_HH
#define G4ATLASRK4_HH

#include "globals.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"

class G4NystromRK4 : public G4MagIntegratorStepper
{
  public: 
    G4NystromRK4(G4EquationOfMotion *EquationMotion, G4int numberOfVariables = 6);
    // G4NystromRK4(G4Mag_EqRhs *EquationMotion);  // Enforces MagEq, numVar=6

    ~G4NystromRK4() ;

    void Stepper(const G4double P   [],
	         const G4double dPdS[],
	               G4double step  ,
	               G4double Po  [],
	               G4double Err []);

    G4int     IntegratorOrder() const {return 4;}
    G4double  DistChord() const; 
    inline void OwnRightHandSide(const double P[],double dPdS[]);   
  
  private:

    void getField   (const G4double P[]);

    ////////////////////////////////////////////////////////////////
    // Private data
    ////////////////////////////////////////////////////////////////

    G4Mag_EqRhs*           m_fEq;          
    G4double      m_field    [3];
    G4double      m_fpos     [4];
    G4double      m_magdistance ;
    G4double      m_magdistance2;
    G4double      m_cof         ;
    G4double      m_mom         ;
    G4double      m_imom        ;
    G4double      m_iPoint   [3];
    G4double      m_mPoint   [3];
    G4double      m_fPoint   [3];
    
};

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// Get new magnetic field with testing distance from previouse call
/////////////////////////////////////////////////////////////////////////////////

inline void G4NystromRK4::getField (const G4double P[])
{
  
  G4double dx = P[0]-m_fpos[0];
  G4double dy = P[1]-m_fpos[1];
  G4double dz = P[2]-m_fpos[2];

  if((dx*dx+dy*dy+dz*dz) > m_magdistance2) {

    m_fpos[0] = P[0];
    m_fpos[1] = P[1];
    m_fpos[2] = P[2];
    m_fpos[3] = P[7];
    m_fEq->GetFieldValue(m_fpos,m_field);
  }
}

/////////////////////////////////////////////////////////////////////////////////
// Derivatives calculation 
/////////////////////////////////////////////////////////////////////////////////

inline void 
G4NystromRK4::OwnRightHandSide(const G4double P[],G4double dPdS[])
{
  getField(P);
  m_mom   = sqrt(P[3]*P[3]+P[4]*P[4]+P[5]*P[5])     ; 
  m_imom  = 1./m_mom                                ;
  m_cof   = m_fEq->FCof()*m_imom                    ;
  dPdS[0] = P[3]*m_imom                             ; // dx /ds
  dPdS[1] = P[4]*m_imom                             ; // dy /ds
  dPdS[2] = P[5]*m_imom                             ; // dz /ds
  dPdS[3] = m_cof*(P[4]*m_field[2]-P[5]*m_field[1]) ; // dPx/ds
  dPdS[4] = m_cof*(P[5]*m_field[0]-P[3]*m_field[2]) ; // dPy/ds
  dPdS[5] = m_cof*(P[3]*m_field[1]-P[4]*m_field[0]) ; // dPz/ds
}


#endif  // G4ATLASRK4_HH
