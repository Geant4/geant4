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
// G4BorisScheme implementation
//
// Author: Divyansh Tiwari, Google Summer of Code 2022
// Supervision: John Apostolakis,Renee Fatemi, Soon Yung Jun
// --------------------------------------------------------------------

#include "G4BorisScheme.hh"
#include "G4FieldUtils.hh"
#include"G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"

#include "G4EquationOfMotion.hh"
//#include "G4EqMagElectricField.hh"

using namespace field_utils;

G4BorisScheme::G4BorisScheme( G4EquationOfMotion* equation,
                                        G4int nvar )
  : fEquation(equation), fnvar(nvar)
{
  if (nvar <= 0)
  {
    G4Exception("G4BorisScheme::G4BorisScheme()",
                "GeomField0002", FatalException,
                "Invalid number of variables; must be greater than zero!");
  }
}

void G4BorisScheme::DoStep(const G4double restMass,const G4double charge, const G4double yIn[], 
                                 G4double yOut[], G4double hstep) const
{
  G4double yOut1Temp[G4FieldTrack::ncompSVEC];
  G4double yOut2Temp[G4FieldTrack::ncompSVEC];
  
  // Used the scheme described in the following paper:https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/153167/eth-5175-01.pdf?sequence=1
  UpdatePosition(restMass, charge,  yIn, yOut1Temp, hstep/2);
  UpdateVelocity(restMass, charge, yOut1Temp, yOut2Temp, hstep);
  UpdatePosition(restMass, charge, yOut2Temp, yOut, hstep/2);
}

void G4BorisScheme::UpdatePosition(const G4double restMass, const G4double /*charge*/, const G4double yIn[],
                                   G4double yOut[], G4double hstep) const
{
    // Particle information
    copy(yOut, yIn);
    
    // Obtaining velocity
    G4ThreeVector momentum_vec =G4ThreeVector(yIn[3],yIn[4],yIn[5]);
    G4double momentum_mag = momentum_vec.mag();
    G4ThreeVector momentum_dir =(1.0/momentum_mag)*momentum_vec;

    G4double velocity_mag = momentum_mag*(c_l)/(std::sqrt(sqr(momentum_mag) +sqr(restMass)));
    G4ThreeVector velocity = momentum_dir*velocity_mag;

    //Obtaining the time step from the length step

    hstep /= velocity_mag*CLHEP::m;

    // Updating the Position
    for(G4int i = 0; i <3; i++ )
    {
      G4double pos = yIn[i]/CLHEP::m;
      pos += hstep*velocity[i];
      yOut[i] = pos*CLHEP::m;
    }   
}

void G4BorisScheme::UpdateVelocity(const G4double restMass, const G4double charge, const G4double yIn[], 
                                   G4double yOut[], G4double hstep) const
{
   //Particle information
    G4ThreeVector momentum_vec =G4ThreeVector(yIn[3],yIn[4],yIn[5]);
    G4double momentum_mag = momentum_vec.mag();
    G4ThreeVector momentum_dir =(1.0/momentum_mag)*momentum_vec;

    G4double gamma = std::sqrt(sqr(momentum_mag) + sqr(restMass))/restMass;  
    
    G4double mass = (restMass/c_squared)/CLHEP::kg;
    
    //Obtaining velocity
   
    G4double velocity_mag = momentum_mag*(c_l)/(std::sqrt(sqr(momentum_mag) +sqr(restMass)));
    G4ThreeVector velocity = momentum_dir*velocity_mag;

    ////Obtaining the time step from the length step
    
    hstep /= velocity_mag*CLHEP::m; 
       
    // Obtaining the field values
    G4double dydx[G4FieldTrack::ncompSVEC];
    G4double fieldValue[6] ={0,0,0,0,0,0};
    fEquation->EvaluateRhsReturnB(yIn, dydx, fieldValue);
   
    //Initializing Vectors
    G4ThreeVector B;
    G4ThreeVector E;
    copy(yOut, yIn);
    for( G4int i = 0; i < 3; i++)
    {
        E[i] = fieldValue[i+3]/CLHEP::volt*CLHEP::meter;// FIXME - Check Units
        B[i] = fieldValue[i]/CLHEP::tesla;   
    }
    
    //Boris Algorithm
    G4double qd = hstep*(charge/(2*mass*gamma));
    G4ThreeVector h = qd*B;
    G4ThreeVector u = velocity + qd*E;
    G4double h_l = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
    G4ThreeVector s_1 = (2*h)/(1 + h_l);
    G4ThreeVector ud = u + (u + u.cross(h)).cross(s_1);
    G4ThreeVector v_fi = ud +qd*E;
    G4double v_mag = std::sqrt(v_fi.mag2());
    G4ThreeVector v_dir = v_fi/v_mag;
    G4double momen_mag = (restMass*v_mag)/(std::sqrt(c_l*c_l - v_mag*v_mag));
    G4ThreeVector momen = momen_mag*v_dir;

    // Storing the updated momentum
    for(int i = 3; i < 6; i++)
    {
        yOut[i] = momen[i-3];   
    }
}

// ----------------------------------------------------------------------------------

void G4BorisScheme::copy(G4double dst[], const G4double src[]) const
{
  std::memcpy(dst, src, sizeof(G4double) * fnvar);
}

// ----------------------------------------------------------------------------------
// - Methods using the Boris Scheme Stepping to estimate integration error
// ----------------------------------------------------------------------------------
void G4BorisScheme::
StepWithErrorEstimate(const G4double yIn[], G4double restMass, G4double charge, G4double hstep,
                      G4double yOut[], G4double yErr[]) const
{
   // Use two half-steps (comparing to a full step) to obtain output and error estimate
   G4double yMid[G4FieldTrack::ncompSVEC];
   StepWithMidAndErrorEstimate(yIn, restMass, charge, hstep, yMid, yOut, yErr);
}

// ----------------------------------------------------------------------------------

void G4BorisScheme::
StepWithMidAndErrorEstimate(const G4double yIn[],  G4double restMass, G4double charge, G4double hstep,
                                  G4double yMid[], G4double yOut[],   G4double yErr[]
   ) const
{
   G4double halfStep= 0.5*hstep;
   G4double yOutAlt[G4FieldTrack::ncompSVEC];   

   // In a single step
   DoStep(restMass, charge, yIn,  yOutAlt, hstep );

   // Same, and also return mid-point evaluation
   DoStep(restMass, charge, yIn,  yMid, halfStep );
   DoStep(restMass, charge, yMid, yOut, halfStep );

   for( G4int i= 0; i<fnvar; i++ )
   {
      yErr[i] = yOutAlt[i] - yOut[i];
   }
}
