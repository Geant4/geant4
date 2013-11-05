#ifndef G4BiasingTrackData_hh
#define G4BiasingTrackData_hh

class G4Track;
class G4VBiasingOperation;
class G4VBiasingOperator;
class G4BiasingProcessInterface;

class G4BiasingTrackData {
public:
  G4BiasingTrackData(const G4Track* track);
  G4BiasingTrackData(const G4Track*                            track,
		     const G4VBiasingOperation*       birthOperation,
		     const G4VBiasingOperator*         birthOperator,
		     const G4BiasingProcessInterface*   birthProcess);
  ~G4BiasingTrackData();
  
  void SetBirthOperation ( const G4VBiasingOperation*          birthOperation ) { fBirthOperation  = birthOperation;  }
  void SetBirthOperator  ( const G4VBiasingOperator*           birthOperator  ) { fBirthOperator   = birthOperator;   }
  void SetBirthProcess( const G4BiasingProcessInterface* birthProcess) { fBirthProcess = birthProcess; }
  
  const G4Track*                   GetTrack()        const { return fTrack; }
  const G4VBiasingOperation*          GetBirthOperation()  const { return fBirthOperation; }
  const G4VBiasingOperator*           GetBirthOperator()   const { return fBirthOperator; }
  const G4BiasingProcessInterface* GetBirthProcess() const { return fBirthProcess; }
  
private:
  const G4Track*                           fTrack;
  const G4VBiasingOperation*      fBirthOperation;
  const G4VBiasingOperator*        fBirthOperator;
  const G4BiasingProcessInterface*  fBirthProcess;
		     
};

#endif
