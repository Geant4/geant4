#ifndef G4Delete_h
#define G4Delete_h
template <class T>
class G4Delete 
{ 
  public: 
    void operator()(T* aT)
    { 
      if(aT)
      {
        delete aT;
      }	
    }
};

#endif
