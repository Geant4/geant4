// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPrimaryGenerator.hh,v 1.2 1999-11-05 04:16:19 asaim Exp $
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
};

#endif

