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

#ifndef G4UPPEVENT_H
#define G4UPPEVENT_H


class G4UppEvent
{
public:

  G4UppEvent(const G4int PA, const G4int PZ, 
	     const G4int TA, const G4int TZ,
	     const G4double eBeam, const G4double b) : 
    ProjectileA(PA), ProjectileZ(PZ), TargetA(TA), TargetZ(TZ),
    BeamEnergy(eBeam), ImpactParameter(b) {};

  void SetProjectile (const G4int A, const G4int Z) 
    { ProjectileA = A; ProjectileZ = Z; }
  void SetTarget (const G4int A, const G4int Z) 
    { TargetA = A; TargetZ = Z; }
  void SetBeamEnergy (const G4double eBeam) 
    { BeamEnergy = eBeam; }
  
private:

  G4int ProjectileA;
  G4int ProjectileZ;
  G4int TargetA;
  G4int TargetZ;
  G4double BeamEnergy;
  G4double ImpactParameter;

};


#endif // G4UPPEVENT_H
