// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRConst.hh,v 1.4 1999-12-15 14:54:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA, Tue Jul  2 15:30:49 JST 1996
///////////////////////////////////
///// G4FRConst.hh /////
///////////////////////////////////

#if !defined G4_FR_COMMAND_LIST_HH
#define G4_FR_COMMAND_LIST_HH

	//----- Header comment
const   char  FR_G4_PRIM_HEADER        [] = "##G4.PRIM-FORMAT-2.4" ;
//const   char  FR_PHYSICAL_VOLUME_NAME  [] = "#/PhysicalVolumeName" ;
const   char  FR_PHYSICAL_VOLUME_NAME  [] = "#/PVName" ;

	//----- Execution control (beginning with !)
const   char  FR_GUI                   [] = "!GraphicalUserInterface" ;
const   char  FR_DEVICE                [] = "!Device"           ;
const   char  FR_SET_CAMERA            [] = "!SetCamera"        ;
const	char  FR_OPEN_DEVICE           [] = "!OpenDevice"       ;
const	char  FR_CLOSE_DEVICE          [] = "!CloseDevice"      ;
const   char  FR_DRAW_ALL              [] = "!DrawAll"          ;
const   char  FR_CLEAR_DATA            [] = "!ClearData"        ;
const   char  FR_QUIT                  [] = "!Quit"             ;
const	char  FR_DISCONNECT_DAWND      [] = "!DisconnectDawnd"  ;
const	char  FR_TERMINATE_DAWND       [] = "!TerminateDawnd"   ;
const   char  FR_SAVE                  [] = "!Save"             ;
const   char  FR_END_SAVE              [] = "!EndSave"          ;
const   char  FR_WAIT                  [] = "!Wait"             ;
const   char  FR_PAUSE                 [] = "!Pause"            ;

	//----- Drawing Style
const	char  FR_WIREFRAME   [] = "/Wireframe" ;
const	char  FR_SURFACE     [] = "/Surface"   ;
const   char  FR_LINES       [] = "/Lines"     ;

	//----- Begin and End of Modeling
const	char  FR_BEGIN_MODELING        [] = "!BeginModeling"    ;
const	char  FR_END_MODELING          [] = "!EndModeling"      ;

	//----- Bounding box
const	char  FR_BOUNDING_BOX   [] = "/BoundingBox"    ;
const	char  FR_BOUNDING_BOX_UNIT [] = "/BoundingBox -0.5 -0.5 -0.5  0.5 0.5 0.5" ;

	//----- Number of divising curved surface
const	int   FR_DEFALUT_NDIV_VALUE = 24 ;

const	char  FR_NDIV           [] = "/Ndiv"           ;
const	char  FR_NDIV_DEFAULT [] = "/Ndiv  24" ;
const	char  FR_NDIV_3  [] = "/Ndiv   3"      ;
const	char  FR_NDIV_4  [] = "/Ndiv   4"      ;
const	char  FR_NDIV_8  [] = "/Ndiv   8"      ;
const	char  FR_NDIV_16 [] = "/Ndiv  16"      ;
const	char  FR_NDIV_24 [] = "/Ndiv  24"      ;
const	char  FR_NDIV_48 [] = "/Ndiv  48"      ;
const	char  FR_NDIV_96 [] = "/Ndiv  96"      ;

	//----- Camera information
const	char  FR_CAMERA_POSITION[] = "/CameraPosition" ;
const	char  FR_CAMERA_POSITION_DEFAULT[] = "/CameraPosition  100000.0  0.0  0.0" ;
				// see from far upward position
const	char  FR_TARGET_POINT   [] = "/TargetPoint";
const	char  FR_ZOOM_FACTOR    [] = "/ZoomFactor";
const	char  FR_SCALE          [] = "/Scale";
const	char  FR_FOCAL_DISTANCE [] = "/FocalDistance";

	//----- Body coordinate information
const	char  FR_BASE_VECTOR    [] = "/BaseVector" ; 
			// Give e1 and e2. Then e3 is calculated. 
const	char  FR_BASE_VECTOR_DEFAULT    [] = "/BaseVector 1.0 0.0 0.0  0.0 1.0 0.0"  ;
const	char  FR_ORIGIN         [] = "/Origin"         ;
const	char  FR_ORIGIN_DEFAULT [] = "/Origin  0.0  0.0  0.0"  ;

	//----- Attribute information
const	char  FR_DIFFUSE_RGB       []  = "/DiffuseRGB" ;
					// old name of /ColorRGB

const	char  FR_COLOR_RGB         []  = "/ColorRGB"   ;
const	char  FR_COLOR_RGB_DEFAULT []  = "/ColorRGB  1.0  1.0  1.0" ;
const	char  FR_COLOR_RGB_WHITE   []  = "/ColorRGB  1.0  1.0  1.0" ;

const	char  FR_COLOR_RGB_RED     []  = "/ColorRGB  1.0  0.0  0.0" ;
const	char  FR_COLOR_RGB_GREEN   []  = "/ColorRGB  0.0  1.0  0.0" ;
const	char  FR_COLOR_RGB_BLUE    []  = "/ColorRGB  0.0  0.0  1.0" ;

const	char  FR_COLOR_RGB_CYAN    []  = "/ColorRGB  0.0  1.0  1.0" ;
const	char  FR_COLOR_RGB_MAGENTA []  = "/ColorRGB  1.0  0.0  1.0" ;
const	char  FR_COLOR_RGB_YELLOW  []  = "/ColorRGB  1.0  1.0  0.0" ;

const	char  FR_SPECULAR_RGB   [] = "/SpecularRGB"  ;
const	char  FR_SPECULAR_RGB_DEFAULT [] = "/SpecularRGB  1.0  1.0  1.0"  ;
const	char  FR_SPECULAR_RGB_WHITE   [] = "/SpecularRGB  1.0  1.0  1.0"  ;

const	char  FR_PHONG_POWER    [] = "/PhongPower"   ;
const	char  FR_PHONG_POWER_DEFAULT [] = "/PhongPower  3"   ;

const	char  FR_TRANSPARENCY   [] = "/Transparency" ;
const	char  FR_TRANSPARENCY_ON   [] = "/Transparency  1 " ; // transparent
const	char  FR_TRANSPARENCY_OFF  [] = "/Transparency  0 " ; // non-transparent

const	char  FR_FORCE_WIREFRAME   [] = "/ForceWireframe" ;
const	char  FR_FORCE_WIREFRAME_ON   [] = "/ForceWireframe  1" ;
const	char  FR_FORCE_WIREFRAME_OFF  [] = "/ForceWireframe  0" ;

const	char  FR_VISIBILITY     [] = "/Visibility"   ;
const	char  FR_VISIBILITY_ON     [] = "/Visibility  1"    ; // visible
const	char  FR_VISIBILITY_OFF    [] = "/Visibility  0"    ; // invisible

	//----- 3D Primitives
const	char  FR_POLYHEDRON     [] = "/Polyhedron"    ;
const	char  FR_VERTEX         [] = "/Vertex"        ;
const	char  FR_FACET          [] = "/Facet"         ;
const	char  FR_END_POLYHEDRON [] = "/EndPolyhedron" ;

const	char  FR_BOX            [] = "/Box"            ; 
const	char  FR_BOX_UNIT  [] = "/Box  0.5  0.5  0.5" ; // dx, dy, dz

const	char  FR_COLUMN         [] = "/Column"         ; 
const	char  FR_COLUMN_UNIT [] = "/Column 0.5  0.5" ; // R dz

const	char  FR_POLYLINE         [] = "/Polyline"          ; 
const	char  FR_PL_VERTEX        [] = "/PLVertex"          ; 
const	char  FR_PL_VERTEX_OLD    [] = "PLVertex"           ; 
const	char  FR_END_POLYLINE     [] = "/EndPolyline"       ; 

const	char  FR_TRD            [] = "/Trd"                  ; 
	// /Trd  dx1 dx2 dy1 dy2 dz    ; 
const	char  FR_TRAP           [] = "/Trap"                 ; 
	// /Trap dz theta phi h1 bl1 tl1 alpha1 h2 bl2 tl2 alpha2
const	char  FR_TUBS           [] = "/Tubs"                 ; 
	// /Tubs rmin rmax dz sphi dphi  
const	char  FR_CONS           [] = "/Cons"                 ; 
	// /Cons rmin1 rmax1 rmin2 rmax2 dz sphi dphi  
const	char  FR_SPHERE           [] = "/Sphere"             ; 
	// /Sphere  R
const	char  FR_SPHERE_SEG       [] = "/SphereSeg"          ; 
	// /SphereSeg  rmin rmax s_theta d_theta s_phi d_phi
const	char  FR_PARA             [] = "/Parallelepiped"     ;
	// /Parallelepiped  dx dy dz tanAlpha tanTheta_cosPhi tanTheta_sinPhi
const	char  FR_PCON             [] = "/PolyCone"           ;
	// /PolyCone  sphi  dphi  nz  z[nz]  rmin[nz]  rmax[nz]
const	char  FR_PGON             [] = "/PolyGon"            ;
	// /PolyGon   sphi  dphi  ndiv  nz  z[nz]  rmin[nz]  rmax[nz]
const	char  FR_TORUS            [] = "/Torus"              ;
	// /PolyGon   sphi  dphi  ndiv  nz  z[nz]  rmin[nz]  rmax[nz]

//----- Marks I (arg: x y z half_size_3d)
const   char  FR_FONT_NAME       [] = "/FontName" ;

//----- Marks I (arg: x y z half_size_3d)
const	char  FR_MARK_CIRCLE_2D  [] = "/MarkCircle2D"   ; 
const	char  FR_MARK_SQUARE_2D  [] = "/MarkSquare2D"   ; 
const	char  FR_MARK_TEXT_2D    [] = "/MarkText2D"   ; 
	// /MarkText2D  x y z  size_world x_offset_world y_offset_world string

//----- Marks II (arg: x y z half_size_2d)
const	char  FR_MARK_CIRCLE_2DS  [] = "/MarkCircle2DS"   ; 
const	char  FR_MARK_SQUARE_2DS  [] = "/MarkSquare2DS"   ; 
const	char  FR_MARK_TEXT_2DS    [] = "/MarkText2DS"   ; 
	// /MarkText2DS  x y z  size_pt x_offset_pt y_offset_pt string

//----- Text
const	char  FR_TEXT_2DS    [] = "/Text2DS"   ; 	

//----- For DAWNCUT
const	char  FR_CLIPPING_PLANE  [] = "/ClippingPlane"   ; 
	// /ClippingPlaneMarkCross2D  a b c d 
	//   for plane   ax + by + cz + d = 0 

#endif
