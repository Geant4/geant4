#ifndef Tst33IStoreBuilder_hh
#define Tst33IStoreBuilder_hh Tst33IStoreBuilder_hh

class G4VIStore;
class Tst33VGeometry;

class Tst33IStoreBuilder {
public:
  Tst33IStoreBuilder();
  ~Tst33IStoreBuilder();
  G4VIStore *CreateIStore(Tst33VGeometry *samplegeo);  
};


#endif
