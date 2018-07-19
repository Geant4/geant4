//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpenGLSceneHandler.hh 99440 2016-09-22 08:34:04Z gcosmo $
//
//
// Andrew Walkden  27th March 1996
// OpenGL scene handler - base for immediate mode and stored mode classes to
//                        inherit from.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSCENEHANDLER_HH
#define G4OPENGLSCENEHANDLER_HH

#include "globals.hh"

#include "G4VSceneHandler.hh"
#include "G4OpenGL.hh"

#include <map>

class G4AttHolder;

// Base class for various OpenGLSceneHandler classes.
class G4OpenGLSceneHandler: public G4VSceneHandler {
  
  friend class G4OpenGLViewer;
  
public:
  virtual void BeginPrimitives (const G4Transform3D& objectTransformation);
  virtual void EndPrimitives ();
  virtual void BeginPrimitives2D (const G4Transform3D& objectTransformation);
  virtual void EndPrimitives2D ();
  void AddPrimitive (const G4Polyline&);
  void AddPrimitive (const G4Polymarker&);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitivesCircle (const std::vector <G4VMarker>&);
  void AddPrimitivesSquare (const std::vector <G4VMarker>&);
  void AddPrimitive (const G4Scale&);
  void AddPrimitive (const G4Polyhedron&);
  
  void PreAddSolid (const G4Transform3D& objectTransformation,
                    const G4VisAttributes&);
  void AddSolid (const G4Box&);
  void AddSolid (const G4Cons&);
  void AddSolid (const G4Tubs&);
  void AddSolid (const G4Trd&);
  void AddSolid (const G4Trap&);
  void AddSolid (const G4Sphere&);
  void AddSolid (const G4Para&);
  void AddSolid (const G4Torus&);
  void AddSolid (const G4Polycone&);
  void AddSolid (const G4Polyhedra&);
  void AddSolid (const G4Orb&);
  void AddSolid (const G4Ellipsoid&);
  void AddSolid (const G4VSolid&);
  void AddCompound (const G4VTrajectory&);
  void AddCompound (const G4VHit&);
  void AddCompound (const G4VDigi&);
  void AddCompound (const G4THitsMap<G4double>&);
  void AddCompound (const G4THitsMap<G4StatDouble>&);
  
  // enum for /vis/ogl/flushAt.
  enum FlushAction {
    endOfEvent,
    endOfRun,
    eachPrimitive,
    NthPrimitive,
    NthEvent,
    never
  };
  static G4int GetEntitiesFlushInterval()
  {return fEntitiesFlushInterval;}
  static FlushAction GetFlushAction()
  {return fFlushAction;}
  static void SetEntitiesFlushInterval(G4int interval)
  {fEntitiesFlushInterval = interval;}
  static void SetFlushAction(FlushAction action)
  {fFlushAction = action;}

#ifdef G4OPENGL_VERSION_2
  private :
  // vertex vector to be given to the graphic card
  std::vector<double> fOglVertex;
  // indices vector to be given to the graphic card
  std::vector<unsigned short> fOglIndices;
  // before, drawyType (as GL_QUADS, GL_TRIANGLES...) was
  // given in glBegin. Now it has to be given in glDrawArray (at the end)
  GLenum fDrawArrayType;
  // emulate GL_QUADS behaviour by inverting two last positions
  bool fEmulate_GL_QUADS;
  // Try to optimize a bit the pipeline
  void OptimizeVBOForTrd();
  void OptimizeVBOForCons(G4int aNoFacet);
  // emulating glEnd and glBegin
  void glEndVBO();
  void glBeginVBO(GLenum type);
  void drawVBOArray(std::vector<double> vertices);
  
  // Buffers used to access vertex and indices elements
#ifndef G4VIS_BUILD_OPENGLWT_DRIVER
  GLuint fVertexBufferObject;
  GLuint fIndicesBufferObject;
#else
  Wt::WGLWidget::Buffer fVertexBufferObject;
  Wt::WGLWidget::Buffer fIndicesBufferObject;
#endif // G4VIS_BUILD_OPENGLWT_DRIVER
  
#endif //G4OPENGL_VERSION_2
  
protected:
  
  G4OpenGLSceneHandler (G4VGraphicsSystem& system,
                        G4int id,
                        const G4String& name = "");
  virtual ~G4OpenGLSceneHandler ();
  
  void ProcessScene();
  G4VSolid* CreateSectionSolid ();
  G4VSolid* CreateCutawaySolid ();
  
  void ClearAndDestroyAtts();  // Destroys att holders and clears pick map.
  
#ifdef G4OPENGL_VERSION_2
  // Special case for VBO, we want to have acces to the VBO drawer everywhere
  // because instead of OpenGL call which are static, VBO openGL functions :
  // - Are functions of an WGLWidget object(G4OpenGLImmediateViewer in our case)
  // - Needs an access to the QGLSHader
  // - Have to be redefined in a VBO way

  inline void setVboDrawer(G4OpenGLVboDrawer* drawer) {
    fVboDrawer = drawer;
  }
  G4OpenGLVboDrawer* fVboDrawer;
#endif // G4OPENGL_VERSION_2
  
  GLuint fPickName;
  std::map<GLuint, G4AttHolder*> fPickMap;  // For picking.
  
  // Shared code to wait until we make a single glFlush
  void ScaledFlush () ;
  // Static so that they apply to all OGL scene handlers...
  static FlushAction fFlushAction;
  // Number of entities between flushes
  static G4int fEntitiesFlushInterval;

  // True if caller of primitives is capable of processing three passes.
  G4bool fThreePassCapable;
  
  G4bool fSecondPassForTransparencyRequested;
  G4bool fSecondPassForTransparency;
  
  G4bool fThirdPassForNonHiddenMarkersRequested;
  G4bool fThirdPassForNonHiddenMarkers;
  
  static const GLubyte fStippleMaskHashed [128];
  bool fEdgeFlag;
};

#include "G4OpenGLSceneHandler.icc"

#endif // G4OPENGLSCENEHANDLER_HH

#endif // G4VIS_BUILD_OPENGL_DRIVER
