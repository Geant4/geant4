#ifndef G4VImportanceScoreConstructor_hh
#define G4VImportanceScoreConstructor_hh G4VImportanceScoreConstructor_hh

class G4VImportanceScoreConstructor {
public:
  virtual ~G4VImportanceScoreConstructor(){}
  virtual void Initialize() = 0;
};

#endif
