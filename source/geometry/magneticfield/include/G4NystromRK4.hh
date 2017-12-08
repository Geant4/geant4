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
// $Id: G4NystromRK4.hh 106564 2017-10-13 09:06:13Z gcosmo $ 
//
// class G4NystromRK4
//
// Class description:
//
// Integrate the equations of the motion of a particle in a magnetic field
// using 4th Runge-Kutta-Nystrom method with errors estimation 
// (ATL-SOFT-PUB-2009-01)
// Current form can be used only for 'pure' magnetic field.
// Notes: 1) field must be time-independent.
//        2) time is not integrated
// 
// History:
// - Created: I.Gavrilenko   15.05.2009   (as G4AtlasRK4)
// - Adaptations:  J. Apostolakis  May-Nov 2009
// -------------------------------------------------------------------

#ifndef G4NYSTROMRK4_HH
#define G4NYSTROMRK4_HH

#include "globals.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"

class G4NystromRK4 : public G4MagIntegratorStepper
{
  public: 
    G4NystromRK4(G4Mag_EqRhs *EquationMotion, G4double distanceConstField=0.0); 
      // Can be used only for Magnetic Fields - and for 6 variables (x,p)

    ~G4NystromRK4() ;

    void Stepper(const G4double P   [],
	         const G4double dPdS[],
	               G4double step  ,
	               G4double Po  [],
	               G4double Err []);
      // Single call for integration result and error
      // - Provides Error via analytical method

    virtual void ComputeRightHandSide(const G4double P[],G4double dPdS[]);   
      // Must compute RHS - and does caches result

    void      SetDistanceForConstantField( G4double length ); 
    G4double  GetDistanceForConstantField() const; 
   
    G4int     IntegratorOrder() const {return 4;}
    G4double  DistChord() const; 
  
  private:

    inline void getField   (const G4double P[4]);

    G4bool CheckCachedMomemtum( const G4double PosMom[6], G4double savedMom );
    G4bool CheckFieldPosition( const G4double Position[3],
                               const G4double lastPosition[3] );
   
    ////////////////////////////////////////////////////////////////
    // Private data
    ////////////////////////////////////////////////////////////////

    G4Mag_EqRhs*           m_fEq;          
    G4double      m_lastField[3];
    G4double      m_fldPosition[4];
    G4double      m_magdistance ;
    G4double      m_magdistance2;
    G4double      m_cof         ;
    G4double      m_mom         ;
    G4double      m_imom        ;
    G4bool        m_cachedMom   ;
    G4double      m_iPoint   [3];
    G4double      m_mPoint   [3];
    G4double      m_fPoint   [3];
    
};

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////
inline void  G4NystromRK4::SetDistanceForConstantField( G4double length )
{
  m_magdistance=   length;
  m_magdistance2 = length*length;
}

inline G4double  G4NystromRK4::GetDistanceForConstantField() const
{
  return m_magdistance; 
}

/////////////////////////////////////////////////////////////////////////////////
// Get value of magnetic field while checking distance from last stored call
/////////////////////////////////////////////////////////////////////////////////

inline void G4NystromRK4::getField (const G4double P[4])
{
  
  G4double dx = P[0]-m_fldPosition[0];
  G4double dy = P[1]-m_fldPosition[1];
  G4double dz = P[2]-m_fldPosition[2];

  if((dx*dx+dy*dy+dz*dz) > m_magdistance2)
  {
    m_fldPosition[0] = P[0];
    m_fldPosition[1] = P[1];
    m_fldPosition[2] = P[2];
    m_fldPosition[3] = P[3];   //  Generally it is P[7] - changed convention !!
    m_fEq->GetFieldValue(m_fldPosition, m_lastField);
  }
}
#endif  // G4NYSTROMRK4
