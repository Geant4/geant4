// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VertexCode.hh,v 1.1.10.1 1999/12/07 20:51:53 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
#ifndef G4VertexCode_h
#define G4VertexCode_h 1


class G4VertexCode
{

  public:
    G4VertexCode(G4String & aCode);

    void SetCode(G4String & aCode);
    G4String GetCode();
  private:

    G4String theCode;
};

inline void G4VertexCode::SetCode(G4String & aCode)
{
  theCode = aCode;
}

inline G4String G4VertexCode::GetCode()
{
  return theCode;
}

#endif
