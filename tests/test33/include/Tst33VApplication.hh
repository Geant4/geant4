#ifndef Tst33VApplication_hh
#define Tst33VApplication_hh Tst33VApplication_hh


class G4VCellScorer;
class G4UserRunAction;
class Tst33VEventAction;

class Tst33VApplication {
public:
  Tst33VApplication();
  virtual ~Tst33VApplication();

  virtual G4UserRunAction *CreateRunAction() = 0;
  virtual Tst33VEventAction *CreateEventAction() = 0;
};

#endif
