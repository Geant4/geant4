#ifndef B03App_hh
#define B03App_hh B03App_hh
class G4RunManager;
class B03PrimaryGeneratorAction;
class B03PhysicsList;
class B03DetectorConstruction;

class B03AppBase {
public:
  ~B03AppBase();
  static B03AppBase &GetB03AppBase();
  G4RunManager &GetRunManager(){return *frunMgr;}
  
private:
  B03AppBase();
  static B03AppBase *fB03AppBase;
  G4RunManager *frunMgr;

  B03DetectorConstruction *fDetector;
  B03PrimaryGeneratorAction *fPrimary;
  B03PhysicsList *fPhysics;
};

// B03AppBase *GetB03AppBase();



#endif










