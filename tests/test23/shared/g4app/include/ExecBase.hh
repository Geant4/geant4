#ifndef ExecBase_H
#define ExecBase_H 1

// fwd declaration
//
class TstReader;
class G4VParticleChange;

class ExecBase
{

   public:
   
      // ctor & dtor
      ExecBase( const TstReader* pset ) { Init(pset); }
      virtual ~ExecBase() {};
      
      
      virtual G4VParticleChange* DoEvent() = 0;

   protected:
   
      ExecBase() {};
      virtual void InitSetup( const TstReader* ) = 0;
      virtual void InitBeam(  const TstReader* ) = 0;
      
      void InitParticles();

   private:
   
      void Init( const TstReader* );

};

#endif
