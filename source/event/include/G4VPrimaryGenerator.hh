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
//
// $Id: G4VPrimaryGenerator.hh,v 1.5 2001-11-20 23:21:41 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VPrimaryGenerator_h
#define G4VPrimaryGenerator_h 1

class G4Event;

// class description:
//
//  This is an abstract base class of all of primary generators.
// This class has only one pure virtual method GeneratePrimaryVertex()
// which takes a G4Event object and generates a primay vertex and
// primary particles associate to the vertex.
//

class G4VPrimaryGenerator
{
  public:
     G4VPrimaryGenerator();
     virtual ~G4VPrimaryGenerator();

     virtual void GeneratePrimaryVertex(G4Event* evt) = 0;

  protected:
     G4ThreeVector         particle_position;
     G4double              particle_time;

  public:
     G4ThreeVector GetParticlePosition()
     { return particle_position; }
     G4double GetParticleTime()
     { return particle_time; }
     void SetParticlePosition(G4ThreeVector aPosition)
     { particle_position = aPosition; }
     void SetParticleTime(G4double aTime)
     { particle_time = aTime; }

};

#endif

