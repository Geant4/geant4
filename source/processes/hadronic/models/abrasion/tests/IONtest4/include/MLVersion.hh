#ifndef MLVersion_h
#define MLVersion_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "RPTofstream.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLVersion
{
  public:
    MLVersion () {
      version    = "v0.5";
      authors    = "Dr F Lei & Dr P Truscott, QinetiQ Ltd, UK";
      coauthors  = "The Geant4 Collaboration";
      sponsor    =
        G4String("European Space Agency, Space Environments and\n") +
        G4String("                                    Effects Analysis Section\n");
#ifdef SPENVIS
      SPENVISver     = "v5.0";
      implementation = "Redhat/Linux Version";
#endif
    };
    ~MLVersion () {};

    G4String version;
    G4String authors;
    G4String coauthors;
    G4String sponsor;

#ifdef SPENVIS
    G4String SPENVISver;
    G4String implementation;
#endif

    friend RPTofstream & operator << (RPTofstream&, MLVersion&);
};
////////////////////////////////////////////////////////////////////////////////
#endif

