// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4InteractionCode.hh,v 1.1 1998/08/22 08:57:15 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef G4InteractionCode_h
#define G4InteractionCode_h 1


class G4InteractionCode
{

  public:
    G4InteractionCode(G4String & aCode);

    void SetCode(G4String & aCode);
    G4String GetCode();

  private:

    G4String theCode;
};

inline void G4InteractionCode::SetCode(G4String & aCode)
{
  theCode = aCode;
}

inline G4String G4InteractionCode::GetCode()
{
  return theCode;
}

#endif
