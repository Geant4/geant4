// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QParentCluster.cc,v 1.1 2000-08-17 13:55:49 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QParentCluster ----------------
//             by Mikhail Kossov, Sept 1999.
//      class fora Parent nuclear cluster in the CHIPS Model
// ----------------------------------------------------------

//#define debug
//#define pdebug

#include "G4QParentClusterVector.hh"

G4QParentCluster::G4QParentCluster(G4int PDGCode): thePDGCode(PDGCode), theProbability(0.){};

G4QParentCluster::G4QParentCluster(G4int PDGCode, G4double prob): 
  thePDGCode(PDGCode), theProbability(prob){};

G4QParentCluster::G4QParentCluster(const G4QParentCluster& rhs)
{
  thePDGCode     = rhs.thePDGCode;
  theProbability = rhs.theProbability;
}

const G4QParentCluster& G4QParentCluster::operator=(const G4QParentCluster& rhs)
{
  thePDGCode     = rhs.thePDGCode;
  theProbability = rhs.theProbability;

  return *this;
}

G4QParentCluster::~G4QParentCluster() {};

// Standard output for G4QParentCluster
ostream& operator<<(ostream& lhs, G4QParentCluster& rhs)
{//      ===============================================
  lhs << "[ParClPDG=" << rhs.GetPDGCode() << ", probab=" << rhs.GetProbability() << "]";
  return lhs;
}

// Standard output for const G4QParentCluster
ostream& operator<<(ostream& lhs, const G4QParentCluster& rhs)
{//      ===============================================
  lhs << "[ParClPDG=" << rhs.GetPDGCode() << ", probab=" << rhs.GetProbability() << "]";
  return lhs;
}





