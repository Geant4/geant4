// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DAWNFILEScene.hh,v 1.1 1999-01-07 16:14:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA


//=================//
#ifdef G4VIS_BUILD_DAWNFILE_DRIVER
//=================//


#ifndef G4DAWNFILE_SCENE_HH
#define G4DAWNFILE_SCENE_HH

#include "globals.hh"

#include "G4VScene.hh"

#include "G4FRofstream.hh"
#include "G4FRConst.hh"


class G4VisAttributes ;
class G4DAWNFILE;


	//-----
class G4DAWNFILEScene: public G4VScene {

public:

	//----- constructor and destructor
  G4DAWNFILEScene (G4DAWNFILE& system, const G4String& name = "");
  ~G4DAWNFILEScene ();

	//----- overriding base class methods
  void AddPrimitive (const G4Polyline& line);
  void AddPrimitive (const G4Polyhedron& p);
  void AddPrimitive (const G4NURBS& nurb);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Polymarker& polymarker) 
       { G4VScene::AddPrimitive (polymarker); }

  virtual void BeginModeling () ; 
  virtual void EndModeling   () { G4VScene::EndModeling ();}

  virtual void BeginPrimitives (const G4Transform3D& objectTransformation);
  virtual void EndPrimitives ();
  void AddThis ( const G4Box&    box    );
  void AddThis ( const G4Cons&   cons   );
  void AddThis ( const G4Tubs&   tubs   );
  void AddThis ( const G4Trd&    trd    );
  void AddThis ( const G4Trap&   trap   );
  void AddThis ( const G4Sphere& sphere );
  void AddThis ( const G4Para&   para   );
  void AddThis ( const G4Torus&  torus  );
  void AddThis ( const G4VSolid& solid  );

  void ClearStore (){}

	//----- public methods inherent to this class
  static G4int GetSceneCount ();
  void         FREndModeling () ;
  G4bool       IsInModeling () { return flag_in_modeling ; }

  G4bool IsSavingG4Prim   ( void ) { return flag_saving_g4_prim ;	}
  void	BeginSavingG4Prim( void ); 
  void	EndSavingG4Prim  ( void ) ;
  void	SetG4PrimFileName() ;

  G4DAWNFILE&  GetSystem   () { return fSystem   ; }
  void         SendBoundingBox   ( void );
  const char*  GetG4PrimFileName () { return fG4PrimFileName ; }


private:

	//----- Utilities etc (common to DAWN and DAWNFILE drivers )
  G4bool    SendVisAttributes ( const G4VisAttributes*  pAV );
  G4bool    IsVisible     ( void ) ;
  G4bool    InitializeFR  ( void ) ;
  void	    SendTransformedCoordinates( void ) ;
  void	    SendPhysVolName           ( void ) ;

	//----- public methods common to DAWN and DAWNFILE drivers
public:

  void	 SendStr   (	const char*	char_string ) ;
  void	 SendStrInt(	const char*	char_string ,
			G4int		ival    );
  void	 SendStrInt3(	const char*	char_string ,
			G4int		ival1  ,
			G4int		ival2  ,
			G4int		ival3   );
  void	 SendStrInt4(	const char*	char_string ,
			G4int		ival1  ,
			G4int		ival2  ,
			G4int		ival3  ,
			G4int		ival4   );
  void   SendStrDouble(	const char*	char_string ,
			G4double	dval   );
  void	 SendStrDouble2(	const char*	char_string ,
				G4double	dval1  ,
				G4double	dval2  );
  void	 SendStrDouble3(	const char*	char_string ,
				G4double	dval1  ,
				G4double	dval2  ,
				G4double	dval3   );

  void	 SendStrDouble4(	const char*	char_string ,
				G4double	dval1  ,
				G4double	dval2  ,
				G4double	dval3  ,
				G4double	dval4  );

  void	 SendStrDouble5(	const char*	char_string ,
				G4double	dval1  ,
				G4double	dval2  ,
				G4double	dval3  ,
				G4double	dval4  ,
				G4double	dval5  );

  void	 SendStrDouble6(	const char*	char_string ,
				G4double	dval1  ,
				G4double	dval2  ,
				G4double	dval3  ,
				G4double	dval4  ,
				G4double	dval5  ,
				G4double	dval6   );

  void   SendStrDouble7(	const char*	char_string ,
				G4double	dval1  ,
				G4double	dval2  ,
				G4double	dval3  ,
				G4double	dval4  ,
				G4double	dval5  ,
				G4double	dval6  ,
				G4double	dval7   );

  void	SendStrDouble11(	const char*	char_string ,
				G4double	dval1  ,
				G4double	dval2  ,
				G4double	dval3  ,
				G4double	dval4  ,
				G4double	dval5  ,
				G4double	dval6  ,
				G4double	dval7  ,
				G4double	dval8  ,
				G4double	dval9  ,
				G4double	dval10  ,
				G4double	dval11   ) ;

  void	 SendIntDouble3(	G4int		ival   ,
				G4double	dval1  ,
				G4double	dval2  ,
				G4double	dval3  );
  void   SendInt3Str(	G4int		ival1  ,
			G4int		ival2  ,
			G4int		ival3  ,
			const char*	char_string );
  void   SendInt4Str(	G4int		ival1  ,
			G4int		ival2  ,
			G4int		ival3  ,
			G4int		ival4  ,
			const char*	char_string );

  void	SendStrDouble6Str(	const char*	char_string1 ,
				G4double	dval1  ,
				G4double	dval2  ,
				G4double	dval3  ,
				G4double	dval4  ,
				G4double	dval5  ,
				G4double	dval6  ,
				const char*	char_string2 );

  void	SendInt   (	G4int 		val );
  void	SendDouble(	G4double 	val );

private:
  G4DAWNFILE&	fSystem;     // Graphics system for this scene.
  static G4int	fSceneIdCount;
  static G4int	fSceneCount;    // No. of existing scenes.

  G4FRofstream	fPrimDest    ;  // defined here
  G4bool	flag_in_modeling ;	
		// true:  FR_BEGIN_MODELING has sent to DAWN, and
		//        FR_END_MODELING   has not sent yet.
		// false:  otherwise
		// 
		// The flag flag_in_modeling is set to "true"
		// in BeginModeling(), and to "false" 
		// in FREndModeling().
		// ( EndModeling() is not used.)

  G4bool	flag_saving_g4_prim ;	

  const int	COMMAND_BUF_SIZE    ;

  char		fG4PrimDestDir [256] ; 
  char          fG4PrimFileName[256] ;
  G4int		fMaxFileNum           ;
	
};

#endif
#endif //G4VIS_BUILD_DAWNFILE_DRIVER
