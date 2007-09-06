#ifndef G4VUSERPHYSICSBIASING_HH
#define G4VUSERPHYSICSBIASING_HH

class G4VUserPhysicsBiasing {

public:
  void Construct();

protected:
  virtual void ConstructBiasing() = 0;
  
private:
  void InitialiseGeneralisedProcessing();

};

inline void G4VUserPhysicsBiasing::Construct()
{
  InitialiseGeneralisedProcessing();
  ConstructBiasing();
}
#endif
