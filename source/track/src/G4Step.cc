// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Step.cc,v 1.1 1999-01-07 16:14:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4Step.cc
//
//  Description:
//    This class represents the Step of a particle tracked.
//    It includes information of 
//      1) List of Step points which compose the Step,
//      2) static information of particle which generated the 
//         Step, 
//      3) trackID and parent particle ID of the Step,
//      4) termination condition of the Step,
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "G4Step.hh"
#include "G4VProcess.hh"

////////////////
G4Step::G4Step()
////////////////
{
  fpPreStepPoint  = new G4StepPoint();
  fpPostStepPoint = new G4StepPoint();
}

/////////////////
G4Step::~G4Step()
/////////////////
{
  delete fpPreStepPoint;
  delete fpPostStepPoint;
}

/////////////////////////////
void G4Step::ShowStep() const
/////////////////////////////
{

// Show header
   G4cout << endl;
   G4cout << "    ++G4Step Information " << endl;
   G4cout.precision(3);

// Show G4Step specific information
   G4cout << "      Address of G4Track    : " << fpTrack << endl;
   G4cout << "      Step Length (mm)      : " << fpTrack->GetStepLength() << endl;
   G4cout << "      Energy Deposit (MeV)  : " << fTotalEnergyDeposit << endl;

// Show G4StepPoint specific information
   G4cout << "      -------------------------------------------------------" 
        << "----------------" <<  endl;
   G4cout << "        StepPoint Information  " << setw(20) << "PreStep" 
                                             << setw(20) << "PostStep" << endl;
   G4cout << "      -------------------------------------------------------" 
        << "----------------" <<  endl;
   G4cout << "         Position - x (mm)   : " 
        << setw(20) << fpPreStepPoint->GetPosition().x() 
        << setw(20) << fpPostStepPoint->GetPosition().x() << endl;
   G4cout << "         Position - y (mm)   : " 
        << setw(20) << fpPreStepPoint->GetPosition().y() 
        << setw(20) << fpPostStepPoint->GetPosition().y() << endl;
   G4cout << "         Position - z (mm)   : " 
        << setw(20) << fpPreStepPoint->GetPosition().z() 
        << setw(20) << fpPostStepPoint->GetPosition().z() << endl;
   G4cout << "         Global Time (ns)    : " 
        << setw(20) << fpPreStepPoint->GetGlobalTime()
        << setw(20) << fpPostStepPoint->GetGlobalTime() << endl;
   G4cout << "         Local Time (ns)     : " 
        << setw(20) << fpPreStepPoint->GetLocalTime() 
        << setw(20) << fpPostStepPoint->GetLocalTime() << endl;
   G4cout << "         Proper Time (ns)    : " 
        << setw(20) << fpPreStepPoint->GetProperTime()
        << setw(20) << fpPostStepPoint->GetProperTime() << endl;
   G4cout << "         Momentum Direct - x : " 
        << setw(20) << fpPreStepPoint->GetMomentumDirection().x()
        << setw(20) << fpPostStepPoint->GetMomentumDirection().x() << endl;
   G4cout << "         Momentum Direct - y : " 
        << setw(20) << fpPreStepPoint->GetMomentumDirection().y()
        << setw(20) << fpPostStepPoint->GetMomentumDirection().y() << endl;
   G4cout << "         Momentum Direct - z : " 
        << setw(20) << fpPreStepPoint->GetMomentumDirection().z()
        << setw(20) << fpPostStepPoint->GetMomentumDirection().z() << endl;
   G4cout << "         Momentum - x (MeV/c): " 
        << setw(20) << fpPreStepPoint->GetMomentum().x()
        << setw(20) << fpPostStepPoint->GetMomentum().x() << endl;
   G4cout << "         Momentum - y (MeV/c): " 
        << setw(20) << fpPreStepPoint->GetMomentum().y()
        << setw(20) << fpPostStepPoint->GetMomentum().y() << endl;
   G4cout << "         Momentum - z (MeV/c): " 
        << setw(20) << fpPreStepPoint->GetMomentum().z()
        << setw(20) << fpPostStepPoint->GetMomentum().z() << endl;
   G4cout << "         Total Energy (MeV)  : " 
        << setw(20) << fpPreStepPoint->GetTotalEnergy()
        << setw(20) << fpPostStepPoint->GetTotalEnergy() << endl;
   G4cout << "         Kinetic Energy (MeV): " 
        << setw(20) << fpPreStepPoint->GetKineticEnergy()
        << setw(20) << fpPostStepPoint->GetKineticEnergy() << endl;
   G4cout << "         Velocity (mm/ns)    : " 
        << setw(20) << fpPreStepPoint->GetVelocity()
        << setw(20) << fpPostStepPoint->GetVelocity() << endl;
   G4cout << "         Volume Name         : " 
        << setw(20) << fpPreStepPoint->GetPhysicalVolume()->GetName()
        << setw(20) << fpPostStepPoint->GetPhysicalVolume()->GetName() << endl;
   G4cout << "         Safety (mm)         : " 
        << setw(20) << fpPreStepPoint->GetSafety()
        << setw(20) << fpPostStepPoint->GetSafety() << endl;
   G4cout << "         Polarization - x    : " 
        << setw(20) << fpPreStepPoint->GetPolarization().x()
        << setw(20) << fpPostStepPoint->GetPolarization().x() << endl;
   G4cout << "         Polarization - y    : " 
        << setw(20) << fpPreStepPoint->GetPolarization().y()
        << setw(20) << fpPostStepPoint->GetPolarization().y() << endl;
   G4cout << "         Polarization - Z    : " 
        << setw(20) << fpPreStepPoint->GetPolarization().z()
        << setw(20) << fpPostStepPoint->GetPolarization().z() << endl;
   G4cout << "         Weight              : " 
        << setw(20) << fpPreStepPoint->GetWeight()
        << setw(20) << fpPostStepPoint->GetWeight() << endl;
   G4cout << "         Step Status         : " ;
        G4StepStatus  tStepStatus = fpPreStepPoint->GetStepStatus();
        if( tStepStatus == fGeomBoundary ){
           G4cout << setw(20) << "Geom Limit";
        } else if ( tStepStatus == fAlongStepDoItProc ){
          G4cout << setw(20) << "AlongStep Proc.";
        } else if ( tStepStatus == fPostStepDoItProc ){
           G4cout << setw(20) << "PostStep Proc";
        } else if ( tStepStatus == fAtRestDoItProc ){
           G4cout << setw(20) << "AtRest Proc";
        } else if ( tStepStatus == fUndefined ){
           G4cout << setw(20) << "Undefined";
        }

        tStepStatus = fpPostStepPoint->GetStepStatus();
        if( tStepStatus == fGeomBoundary ){
           G4cout << setw(20) << "Geom Limit";
        } else if ( tStepStatus == fAlongStepDoItProc ){
           G4cout << setw(20) << "AlongStep Proc.";
        } else if ( tStepStatus == fPostStepDoItProc ){
           G4cout << setw(20) << "PostStep Proc";
        } else if ( tStepStatus == fAtRestDoItProc ){
           G4cout << setw(20) << "AtRest Proc";
        } else if ( tStepStatus == fUndefined ){
           G4cout << setw(20) << "Undefined";
        }

        G4cout << endl;
        G4cout << "         Process defined Step: " ;
        if( fpPreStepPoint->GetProcessDefinedStep() == NULL ){
 	   G4cout << setw(20) << "Undefined";
        } else {
  	   G4cout << setw(20) << fpPreStepPoint->GetProcessDefinedStep()
                                             ->GetProcessName();
        }
        if( fpPostStepPoint->GetProcessDefinedStep() == NULL){
  	   G4cout << setw(20) << "Undefined";
        } else {
 	   G4cout << setw(20) << fpPostStepPoint->GetProcessDefinedStep()
                                              ->GetProcessName(); 
        }

   G4cout << endl;
   G4cout << "      -------------------------------------------------------" 
        << "----------------" << endl;
}





