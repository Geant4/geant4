#ifndef G4VPScorer_hh
#define G4VPScorer_hh G4VPScorer_hh

  
class G4Step;
class G4PStep;
  
class G4VPScorer {
public:
  virtual ~G4VPScorer(){}
  virtual void Score(const G4Step &step, const G4PStep &pstep) = 0;
  
};
  

#endif
