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
// G4BorisSDC implementation
//
// Author: Divyansh Tiwari, Google Summer of Code 2022
// Supervision: John Apostolakis,Renee Fatemi, Soon Yung Jun
// --------------------------------------------------------------------

#include "G4BorisSDC.hh"
#include "G4FieldUtils.hh"
#include"G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
using namespace field_utils;

G4BorisSDC::G4BorisSDC( G4EquationOfMotion* equation,
                                        G4int nvar )
  : fEquation(equation), fnvar(nvar)
{
  if (nvar <= 0)
  {
    G4Exception("G4BorisSDC::G4BorisSDC()",
                "GeomField0002", FatalException,
                "Invalid number of variables; must be greater than zero!");
  }
}



G4ThreeVector G4BorisSDC::GetLorentzForce(const G4ThreeVector position, const G4ThreeVector velocity) const
{
    // initialize variables
    G4double dydx[G4FieldTrack::ncompSVEC];
    G4double y[G4FieldTrack::ncompSVEC];

    //convert velocity into momentum
    G4double v_mag = velocity.mag();
    G4ThreeVector v_dir = velocity/v_mag;
    G4double momen_mag = (restMass_c2*v_mag)/(std::sqrt(c_l*c_l - v_mag*v_mag));
    G4ThreeVector momen = momen_mag*v_dir;

    // obtain the y-array
    for(int i = 0; i <3; i++)
    {
        y[i] = position[i];
        y[i+2] = momen[i];
    }
    
    //Calculate field
    G4double fieldValue[6] ={0,0,0,0,0,0};
    fEquation->EvaluateRhsReturnB(y, dydx, fieldValue);

    //initialize E and B
    G4ThreeVector B;
    G4ThreeVector E;
    for( G4int i = 0; i < 3; i++)
    {
        E[i] = fieldValue[i+3]/CLHEP::volt*CLHEP::meter;// FIXME - Check Units
        B[i] = fieldValue[i]/CLHEP::tesla;   
    }

    G4ThreeVector force = alpha*(E + velocity.cross(B));

   // G4cout<<force[0]<<" "<<force[1]<<" "<<force[2]<<G4endl;
    return force; 
}

void G4BorisSDC::DoStep(const G4double restMass,const G4double charge, const G4double yIn[], 
                                 G4double yOut[], G4double hstep) 
{
    // Used the scheme described in the following paper:https://arxiv.org/pdf/1409.5677.pdf

    //Initialize mass and charge
   
    mass_si = (restMass/c_squared)/CLHEP::kg;
    restMass_c2 = restMass;
    G4cout<<"YIN: "<<yIn[0]<<" "<<yIn[1]<<" "<<yIn[2]<<" "<<yIn[3]<<" "<<yIn[4]<<" "<<yIn[5]<<G4endl;
    charge_si = charge;
    alpha = charge_si/mass_si;
   // G4cout<<"Alpha: "<<alpha<<G4endl;

    //initial velocity and position
    G4ThreeVector momentum_vec = G4ThreeVector(yIn[3],yIn[4],yIn[5]);
    G4double momentum_mag = momentum_vec.mag();
    G4ThreeVector momentum_dir =(1.0/momentum_mag)*momentum_vec;

    G4double velocity_mag = momentum_mag*(c_l)/(std::sqrt(sqr(momentum_mag) +sqr(restMass_c2)));
    G4ThreeVector velocity = momentum_dir*velocity_mag;
    G4ThreeVector position = G4ThreeVector(yIn[0], yIn[1], yIn[2]);
    //G4cout<<"Lorentz Force "<<GetLorentzForce(position, velocity)<<G4endl;

    //Get Time step
    hstep /= velocity_mag;

    // calculate delta_t_m a.k.a difference between nodes
    
    //G4cout<<"Delta tm"<<G4endl;
    for(int i = 1; i <= M; i++)
    {
        delta_t_m[i] =hstep*(nodes[i] - nodes[i-1]);
        //G4cout<<delta_t_m[i]<<" ";
    }
   // G4cout<<G4endl;
   // G4cout<<G4endl;

   

    //Calculate S Matrix
    for(int i = 1; i<=M; i++)
    {
        for(int j = 1; j<=M; j++)
        {
           
            if (i == 1) S[i][j] = Q[i][j];
            else  S[i][j] = Q[i][j] - Q[i-1][j];

            S[i][j] *= hstep;

            //G4cout<<S[i][j]<<" ";
        }
        //G4cout<<G4endl;
    }

    // Initialise Velocity and Position matrices
    //G4cout<<"Initial Values"<<G4endl;
    for(int i = 0; i <=M; i++)
    {
        Velocity[i][0] = velocity;
        Position[i][0] = position;

        if(i == 0)
        {
            for( int j = 1; j <=K; j++)
            {
                Velocity[i][j] = velocity;
                Position[i][j] = position;

                //G4cout<<Velocity[i][i][0]<<Velocity[i][i][1]<<Velocity[i][i][2]<<Position[i][j][0]<<Position[i][j][1]<<Position[i][j][2]<<G4endl;
            }
        }
         //G4cout<<"Velocity "<< velocity<<" Position "<<position<<G4endl;


    }

    G4cout<<"Position and Velocity Updates"<<G4endl;
    for(int i = 1; i<=M; i++)
    {
        for(int j = 1; j <=K; j++)
        {

            //Update Position
            UpdatePosition(j,i);
            //Update Velocity
            UpdateVelocity(j,i);
            

           G4cout<<Position[i][j][0]<<" "<<Position[i][j][1]<<" "<<Position[i][j][2]<<" "<<Velocity[i][j][0]<<" "<<Velocity[i][j][1]<<" "<<Velocity[i][j][2]<<G4endl;
        }
       
    }
    G4cout<<G4endl;


    // update the values
    copy(yOut, yIn);
    G4ThreeVector v = Velocity[M][K];
    //convert velocity into momentum
    G4double v_mag = v.mag();
    G4ThreeVector v_dir = v/v_mag;
    G4double momen_mag = (restMass_c2*v_mag)/(std::sqrt(c_l*c_l - v_mag*v_mag));
    G4ThreeVector momen = momen_mag*v_dir;

    for(int i = 0; i < 3; i++)
    {
        yOut[i] = Position[M][K][i];
        yOut[i+3] = momen[i];

    }

    G4cout<<"YOUT: "<<yOut[0]<<" "<<yOut[1]<<" "<<yOut[2]<<" p:  "<<yOut[3]<<" "<<yOut[4]<<" "<<yOut[5]<<G4endl<<G4endl;
}

void G4BorisSDC::UpdatePosition(G4int k, G4int m_i) 
{
    G4ThreeVector param1 = S[m_i][1]*Velocity[1][k-1]; //  IV term in the derivation
    for(int j = 2; j <=M; j++)
    {
        param1 += S[m_i][j]*Velocity[j][k-1];
    }

   G4ThreeVector param2 = Velocity[m_i-1][k] + delta_t_m[m_i]*GetLorentzForce(Position[m_i-1][k], Velocity[m_i-1][k])/2; //half-node velocity 1 
   G4ThreeVector param3 = Velocity[m_i-1][k-1] + delta_t_m[m_i]*GetLorentzForce(Position[m_i-1][k-1], Velocity[m_i-1][k-1])/2; //half-node velocity 2

   Position[m_i][k] = Position[m_i-1][k] + delta_t_m[m_i]*(param2 - param3) + param1; // calculate position for this step  
}

void G4BorisSDC::UpdateVelocity( G4int k, G4int m_i) 
{
    G4ThreeVector c_k = -(GetLorentzForce(Position[m_i][k-1], Velocity[m_i][k-1]) + GetLorentzForce(Position[m_i-1][k-1], Velocity[m_i-1][k-1]))/2; // add magnetic field terms later for non-constant magnetic field (Eq 53-54) : alpha/2 *(Vn x (Bn - Bn-1))
    //c_k /= delta_t_m[m_i];
   // G4cout<<"c_k Before : "<<c_k[0]<<" "<<c_k[1]<<" "<<c_k[2]<<G4endl;
   //  G4cout<<"delta_t"<<delta_t_m[m_i]<<G4endl;
    //c_k *= delta_t_m[m_i];
    for(int l =  1; l <= M; l++)
    {
        G4ThreeVector temp =  S[m_i][l]*GetLorentzForce(Position[l][k-1], Velocity[l][k-1]);
        c_k += temp;
       // G4cout<<"c_k iteration:"<<l<<": "<<temp<<" S: "<<S[m_i][l];
    }
   // G4cout<<"c_k  : "<<c_k[0]<<" "<<c_k[1]<<" "<<c_k[2]<<G4endl;
    //G4cout<<G4endl;

   // c_k /= delta_t_m[m_i];


    // Calculate E and B
    G4ThreeVector E,B;
    G4double dydx[G4FieldTrack::ncompSVEC];
    G4double y[G4FieldTrack::ncompSVEC];

    for(int i = 0; i <3; i++)
    {
        y[i] = Position[m_i][k][i];
        y[i+2] = 0;
    }

    //Calculate field
    G4double fieldValue[6] ={0,0,0,0,0,0};
    fEquation->EvaluateRhsReturnB(y, dydx, fieldValue);

    //initialize E and B
    for( G4int i = 0; i < 3; i++)
    {
        E[i] = fieldValue[i+3]/CLHEP::volt*CLHEP::meter;// FIXME - Check Units
        B[i] = fieldValue[i]/CLHEP::tesla;   
    }


    // Calculating Velocity via the rotation scheme defined for a normal boris algorithm
    G4ThreeVector v_minus = Velocity[m_i-1][k] + delta_t_m[m_i]*(alpha*E + c_k)/2;
    G4double gamma; // how to find gamma
    G4ThreeVector t = alpha*B*delta_t_m[m_i]/2;
    G4double t_l = t[0]*t[0] + t[1]*t[1]+ t[2]*t[2];
    G4ThreeVector s_i = 2*t/(1 + t_l);
    G4ThreeVector v_plus = v_minus + (v_minus + v_minus.cross(t)).cross(s_i);
    G4ThreeVector v = v_plus + delta_t_m[m_i]*(alpha*E + c_k)/2;
    //G4cout<<"c_k : "<<c_k[0]<<" "<<c_k[1]<<" "<<c_k[2]<<G4endl;
   // G4cout<<"V_minus: "<<v_minus[0]<<" "<<v_minus[1]<<" v_plus: "<<v_plus[0]<<" "<<v_plus[1]<<G4endl;

    Velocity[m_i][k] = v; 
}


void G4BorisSDC::copy(G4double dst[], const G4double src[]) const
{
  std::memcpy(dst, src, sizeof(G4double) * fnvar);
}
