#ifndef PCTWriter_hh
#define PCTWriter_hh 1

#include <iostream>
#include <iomanip>

#include "gzstream.h"

#include "globals.hh"
#include "G4ReactionProductVector.hh"

class PCTCompositeNucleus;
class G4Fragment;

class PCTWriter
{
public:
  inline PCTWriter();
  inline ~PCTWriter();

  inline void OpenFile(const G4String name);
  inline G4bool IsOpen() const;
  inline void CloseFile();

  inline void WriteString(const G4String& str);

  void WriteHeader(const PCTCompositeNucleus& aCN, const G4double aXS);
  void WriteReaction(const G4int iteration, const G4ReactionProductVector * theDynamicPartVectorPtr,
		     const G4Fragment * theExcitedNucleus);

private:
  inline PCTWriter(const PCTWriter& right);
  inline PCTWriter& operator=(const PCTWriter& right);

private:
  ogzstream theFile; 

};

#include "PCTWriter.icc"

#endif // PCTWriter_hh
