#ifndef B01VApplication_hh
#define B01VApplication_hh B01VApplication_hh

class B01VSimulation;

class B01VApplication {
public:
  B01VApplication();
  virtual ~B01VApplication();
  virtual void RunSimulation(B01VSimulation *sim) = 0;
};

#endif
