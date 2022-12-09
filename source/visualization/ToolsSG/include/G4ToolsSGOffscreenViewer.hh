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
// Guy Barrand 09th June 2022

#ifndef G4TOOLSSGOFFSCREENVIEWER_HH
#define G4TOOLSSGOFFSCREENVIEWER_HH

#include "G4ToolsSGViewer.hh"

#include "G4UIcmdWithABool.hh"

#include <tools/fpng>
#include <tools/toojpeg>

#include <tools/offscreen/sg_viewer>
#include <tools/tos>
#include <tools/sto>

class G4ToolsSGOffscreenViewer : public G4ToolsSGViewer<tools::offscreen::session,tools::offscreen::sg_viewer> {
  typedef G4ToolsSGViewer<tools::offscreen::session,tools::offscreen::sg_viewer> parent;
public:
  G4ToolsSGOffscreenViewer(tools::offscreen::session& a_session,G4ToolsSGSceneHandler& a_scene_handler, const G4String& a_name)
  :parent(a_session,a_scene_handler,a_name)
  ,fFileName("auto")
  ,fFilePrefix("auto")
  ,fFileIndex(0)
  ,fResetFileIndex(false)
  {
    Messenger::Create();
  }
  virtual ~G4ToolsSGOffscreenViewer() = default;
protected:
  G4ToolsSGOffscreenViewer(const G4ToolsSGOffscreenViewer& a_from):parent(a_from){}
  G4ToolsSGOffscreenViewer& operator=(const G4ToolsSGOffscreenViewer&) {return *this;}
public:  
  virtual void Initialise() {
    if(fSGViewer) return; //done.
    //::printf("debug : G4ToolsSGOffscreenViewer::Initialize\n");
    // auto refresh true produces too much files.
    fVP.SetAutoRefresh(false);
    fDefaultVP.SetAutoRefresh(false);
    fSGViewer = new tools::offscreen::sg_viewer(fSGSession
      ,fVP.GetWindowAbsoluteLocationHintX(1440)
      ,fVP.GetWindowAbsoluteLocationHintY(900)
      ,fVP.GetWindowSizeHintX()
      ,fVP.GetWindowSizeHintY()
      ,fName);
    fSGViewer->set_file_format("zb_png");
    fSGViewer->set_file_name("out.png");
    fSGViewer->set_png_writer(tools::fpng::write);
    fSGViewer->set_jpeg_writer(tools::toojpeg::write);
    fSGViewer->set_do_transparency(true);
    fSGViewer->set_top_to_bottom(false); //if using tools::fpng, tools::toojpeg.
  }    
  virtual void SetView() {
    //::printf("debug : G4ToolsSGOffscreenViewer::SetView\n");
    fVP.SetGlobalMarkerScale(1);  //WARNING: for __APPLE__, the G4ToolsSGQtViewer set it to 2.
    parent::SetView();
  }

  virtual void DrawView() {
    if (!fNeedKernelVisit) KernelVisitDecision();
    fLastVP = fVP;
    ProcessView();  // Clears store and processes scene only if necessary.
    //::printf("debug : G4ToolsSGOffscreenViewer::DrawView %s\n",fName.c_str());
    if(fSGViewer) {
      fSGSceneHandler.TouchPlotters(fSGViewer->sg());
      if(fFileName=="auto") {
        std::string prefix;
        if(fFilePrefix=="auto") {
          prefix = "g4tsg_offscreen_"+fSGViewer->file_format()+"_";
        } else {
          prefix = fFilePrefix;
        }
        std::string suffix;
        if(G4ToolsSGOffscreenViewer::GetFormatExtension(fSGViewer->file_format(),suffix)) {
          std::string file_name = prefix+tools::tos(GetFileIndex())+"."+suffix;
          fSGViewer->set_file_name(file_name);
        }
      } else {
        fSGViewer->set_file_name(fFileName);
      }
      if(fSGViewer->write_paper()) {
        if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations) {
          G4cout << "File " << fSGViewer->file_name() << " produced." << G4endl;
        }
      }
    }
  }

  virtual void ClearView() {
    //::printf("debug : G4ToolsSGOffscreenViewer::ClearView %s\n",fName.c_str());
  }
  virtual void ShowView() {
    //::printf("debug : G4ToolsSGOffscreenViewer::ShowView %s\n",fName.c_str());
  }
  virtual void FinishView() {
    //::printf("debug : G4ToolsSGOffscreenViewer::FinishView %s\n",fName.c_str());
    if(fSGViewer) {
      fSGSceneHandler.TouchPlotters(fSGViewer->sg());
    }
  }

protected:  
  void SetSize(unsigned int a_w,unsigned int a_h) {
    if(!fSGViewer) return;
    if(!a_w || !a_h) {
      fSGViewer->set_size(fVP.GetWindowSizeHintX(),fVP.GetWindowSizeHintY());
    } else {
      fSGViewer->set_size(a_w,a_h);
    }
  }

  void SetFileFormat(const G4String& a_format) {
    if(!fSGViewer) return;
    fSGViewer->set_file_format(a_format);
  } 
    
  void SetDoTransparency(bool a_value) {
    if(!fSGViewer) return;
    fSGViewer->set_do_transparency(a_value);
  }

  void SetFileName(const G4String& a_file,const G4String& a_prefix,bool a_reset_index) {
    fFileName = a_file;
    fFilePrefix = a_prefix;
    fResetFileIndex = a_reset_index;
  } 

  void SetGL2PSSort(const G4String& a_sort) {
    if(!fSGViewer) return;
    fSGViewer->set_opts_1(a_sort);
  } 
    
  void SetGL2PSOptions(const G4String& a_opts) {
    if(!fSGViewer) return;
    fSGViewer->set_opts_2(a_opts);
  } 
    
  class Messenger: public G4VVisCommand {
  public:  
    static void Create() {static Messenger s_messenger;}
  private:  
    Messenger() {
      G4UIparameter* parameter;
      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      cmd_format = new G4UIcommand("/vis/tsg/offscreen/set/format", this);
      cmd_format->SetGuidance("Set file format.");
      cmd_format->SetGuidance("Available formats are:");
      cmd_format->SetGuidance("- zb_png: tools::sg offscreen zbuffer put in a png file.");
      cmd_format->SetGuidance("- zb_jpeg: tools::sg offscreen zbuffer put in a jpeg file.");
      cmd_format->SetGuidance("- zb_ps: tools::sg offscreen zbuffer put in a PostScript file.");
      cmd_format->SetGuidance("- gl2ps_eps: gl2ps producing eps");
      cmd_format->SetGuidance("- gl2ps_ps:  gl2ps producing ps");
      cmd_format->SetGuidance("- gl2ps_pdf: gl2ps producing pdf");
      cmd_format->SetGuidance("- gl2ps_svg: gl2ps producing svg");
      cmd_format->SetGuidance("- gl2ps_tex: gl2ps producing tex");
      cmd_format->SetGuidance("- gl2ps_pgf: gl2ps producing pgf");

      parameter = new G4UIparameter("format",'s',true);
      parameter->SetDefaultValue("gl2ps_eps");
      cmd_format->SetParameter (parameter);

      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      cmd_file = new G4UIcommand("/vis/tsg/offscreen/set/file", this);
      cmd_file->SetGuidance("Set file name.");
      cmd_file->SetGuidance("Default file name is \"auto\" and default format is zb_png.");
      cmd_file->SetGuidance("If file name is \"auto\", the output file name is built from");
      cmd_file->SetGuidance("a viewer index counter with the form:");
      cmd_file->SetGuidance("    g4tsg_offscreen_<format>_<index>.<format extension>");
      cmd_file->SetGuidance("For example:");
      cmd_file->SetGuidance("    g4tsg_offscreen_zb_png_1.png");
      cmd_file->SetGuidance("    g4tsg_offscreen_zb_png_2.png");
      cmd_file->SetGuidance("    ...");
      cmd_file->SetGuidance("or if format is changed to \"gl2ps_pdf\":");
      cmd_file->SetGuidance("    g4tsg_offscreen_gl2ps_pdf_3.pdf");
      cmd_file->SetGuidance("If a prefix parameter is given, the output file name is built from");
      cmd_file->SetGuidance("a global index counter with the form:");
      cmd_file->SetGuidance("    <prefix><index>.<format extension>");
      cmd_file->SetGuidance("For example:");
      cmd_file->SetGuidance("    /vis/tsg/offscreen/set/file auto my_prefix_");
      cmd_file->SetGuidance("will produce:");
      cmd_file->SetGuidance("    my_prefix_1.png");
      cmd_file->SetGuidance("    my_prefix_2.png");
      cmd_file->SetGuidance("    ...");
      cmd_file->SetGuidance("You can reset the index by specifying true as last argument:");
      cmd_file->SetGuidance("    /vis/tsg/offscreen/set/file auto other_prefix_ true");
      cmd_file->SetGuidance("will produce:");
      cmd_file->SetGuidance("    other_prefix_1.png");
      cmd_file->SetGuidance("    other_prefix_2.png");
      cmd_file->SetGuidance("    ...");

      parameter = new G4UIparameter("file",'s',true);
      parameter->SetDefaultValue("auto");
      cmd_file->SetParameter (parameter);

      parameter = new G4UIparameter("prefix",'s',true);
      parameter->SetDefaultValue("auto");
      cmd_file->SetParameter (parameter);

      parameter = new G4UIparameter("reset_index",'b',true);
      parameter->SetDefaultValue("false");
      cmd_file->SetParameter (parameter);

      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      cmd_size = new G4UIcommand("/vis/tsg/offscreen/set/size", this);
      cmd_size->SetGuidance("Set viewer size in pixels.");
      cmd_size->SetGuidance
	("If width and/or height is set to zero, the viewer size specified with /vis/viewer/create (/vis/open) is taken.");
      cmd_size->SetGuidance(" About the picture size, note that the gl2ps files will grow with the number of primitives");
      cmd_size->SetGuidance("(gl2ps does not have a zbuffer logic). The \"zb\" files will not grow with the number of");
      cmd_size->SetGuidance("primitives, but with the size of the viewer. It should be preferred for scenes with");
      cmd_size->SetGuidance("a lot of objects to render. With zb, to have a better rendering, do not hesitate to");
      cmd_size->SetGuidance("have a large viewer size.");

      parameter = new G4UIparameter("width",'i',false);
      parameter->SetDefaultValue("0");
      cmd_size->SetParameter (parameter);
      
      parameter = new G4UIparameter("height",'i',false);
      parameter->SetDefaultValue("0");
      cmd_size->SetParameter (parameter);

      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      cmd_do_transparency = new G4UIcmdWithABool("/vis/tsg/offscreen/set/transparency", this);
      cmd_do_transparency->SetGuidance("True/false to enable/disable rendering of transparent objects.");
      cmd_do_transparency->SetGuidance("This may be usefull if using file formats, as the gl2ps ones, unable to handle transparency.");
      cmd_do_transparency->SetParameterName("transparency-enabled",true);
      cmd_do_transparency->SetDefaultValue(true);

      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      cmd_gl2ps_sort = new G4UIcommand("/vis/tsg/offscreen/gl2ps/set/sort", this);
      cmd_gl2ps_sort->SetGuidance("Set gl2ps sort algorithm when creating the file.");

      cmd_gl2ps_sort->SetGuidance("The sort argument could be:");
      cmd_gl2ps_sort->SetGuidance(" NO_SORT");
      cmd_gl2ps_sort->SetGuidance(" SIMPLE_SORT");
      cmd_gl2ps_sort->SetGuidance(" BSP_SORT");
      cmd_gl2ps_sort->SetGuidance("The default being BSP_SORT");

      parameter = new G4UIparameter("sort",'s',true);
      parameter->SetDefaultValue("BSP_SORT");
      cmd_gl2ps_sort->SetParameter (parameter);
      
      //////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////
      cmd_gl2ps_opts = new G4UIcommand("/vis/tsg/offscreen/gl2ps/set/options", this);
      cmd_gl2ps_opts->SetGuidance("Set gl2ps options passed when creating the file.");

      cmd_gl2ps_opts->SetGuidance("Options is a list of items separated by |. An item can be:");
      cmd_gl2ps_opts->SetGuidance(" NONE");
      cmd_gl2ps_opts->SetGuidance(" DRAW_BACKGROUND");
      cmd_gl2ps_opts->SetGuidance(" SIMPLE_LINE_OFFSET");
      cmd_gl2ps_opts->SetGuidance(" SILENT");
      cmd_gl2ps_opts->SetGuidance(" BEST_ROOT");
      cmd_gl2ps_opts->SetGuidance(" OCCLUSION_CULL");
      cmd_gl2ps_opts->SetGuidance(" NO_TEXT");
      cmd_gl2ps_opts->SetGuidance(" LANDSCAPE");
      cmd_gl2ps_opts->SetGuidance(" NO_PS3_SHADING");
      cmd_gl2ps_opts->SetGuidance(" NO_PIXMAP");
      cmd_gl2ps_opts->SetGuidance(" USE_CURRENT_VIEWPORT");
      cmd_gl2ps_opts->SetGuidance(" COMPRESS");
      cmd_gl2ps_opts->SetGuidance(" NO_BLENDING");
      cmd_gl2ps_opts->SetGuidance(" TIGHT_BOUNDING_BOX");
      cmd_gl2ps_opts->SetGuidance(" NO_OPENGL_CONTEXT");
      cmd_gl2ps_opts->SetGuidance(" NO_TEX_FONTSIZE");
      cmd_gl2ps_opts->SetGuidance(" PORTABLE_SORT");
      cmd_gl2ps_opts->SetGuidance("The default (typical) list of options is:");
      cmd_gl2ps_opts->SetGuidance(" SILENT|OCCLUSION_CULL|BEST_ROOT|DRAW_BACKGROUND");

      parameter = new G4UIparameter("options",'s',true);
      parameter->SetDefaultValue("SILENT|OCCLUSION_CULL|BEST_ROOT|DRAW_BACKGROUND");
      cmd_gl2ps_opts->SetParameter (parameter);
      

    }
    virtual ~Messenger() {
      delete cmd_format;
      delete cmd_file;
      delete cmd_size;
      delete cmd_do_transparency;
      delete cmd_gl2ps_sort;
      delete cmd_gl2ps_opts;
    }
  public:
    virtual void SetNewValue(G4UIcommand* a_cmd,G4String a_value) {
      G4VisManager::Verbosity verbosity = GetVisManager()->GetVerbosity();
      G4VViewer* viewer = GetVisManager()->GetCurrentViewer();
      if (!viewer) {
        if (verbosity >= G4VisManager::errors) G4cerr << "ERROR: No current viewer." << G4endl;
        return;
      }
      G4ToolsSGOffscreenViewer* tsg_viewer = dynamic_cast<G4ToolsSGOffscreenViewer*>(viewer);
      if(!tsg_viewer) {
        G4cout << "G4ToolsSGOffscreenViewer::Messenger::SetNewValue:"
               << " current viewer is not a G4ToolsSGOffscreenViewer." << G4endl;
        return;
      }
      std::vector<std::string> args;
      tools::double_quotes_tokenize(a_value,args);
      if(args.size()!=a_cmd->GetParameterEntries()) return;
      if(a_cmd==cmd_format) {
	if(!IsKnownFormat(args[0])) {
          G4cout << "G4ToolsSGOffscreenViewer::Messenger::SetNewValue:"
                 << " unknown file format " << args[0] << "." << G4endl;
          return;
	}
        tsg_viewer->SetFileFormat(args[0]);
      } else if(a_cmd==cmd_file) {
        G4bool reset_index = G4UIcommand::ConvertToBool(args[2].c_str());
        tsg_viewer->SetFileName(args[0],args[1],reset_index);
      } else if(a_cmd==cmd_size) {
        unsigned int w,h;
        if(!tools::to(args[0],w)) w = 0;
        if(!tools::to(args[1],h)) h = 0;
        tsg_viewer->SetSize(w,h);
      } else if(a_cmd==cmd_do_transparency) {
        G4bool _do = a_cmd->ConvertToBool(args[0].c_str());
        tsg_viewer->SetDoTransparency(_do);
      } else if(a_cmd==cmd_gl2ps_sort) {
        tsg_viewer->SetGL2PSSort(args[0]);
      } else if(a_cmd==cmd_gl2ps_opts) {
        tsg_viewer->SetGL2PSOptions(args[0]);
      }
    }
  private:
    G4UIcommand* cmd_format;
    G4UIcommand* cmd_file;
    G4UIcommand* cmd_size; 
    G4UIcmdWithABool* cmd_do_transparency;
    G4UIcommand* cmd_gl2ps_sort;
    G4UIcommand* cmd_gl2ps_opts;
 };

protected:  
  unsigned int GetFileIndex() {
    if(fResetFileIndex) {fFileIndex = 0;fResetFileIndex = false;}
    fFileIndex++;
    return fFileIndex;
  }
  
  static bool IsKnownFormat(const std::string& a_format) {
    if(a_format=="gl2ps_eps") return true;
    if(a_format=="gl2ps_ps")  return true;
    if(a_format=="gl2ps_pdf") return true;
    if(a_format=="gl2ps_svg") return true;
    if(a_format=="gl2ps_tex") return true;
    if(a_format=="gl2ps_pgf") return true;
    if(a_format=="zb_ps")     return true;
    if(a_format=="zb_png")    return true;
    if(a_format=="zb_jpeg")   return true;
    return false;
  }

  static bool GetFormatExtension(const std::string& a_format,std::string& a_ext) {
    if(a_format=="gl2ps_eps") {a_ext = "eps";return true;}
    if(a_format=="gl2ps_ps")  {a_ext = "ps";return true;}
    if(a_format=="gl2ps_pdf") {a_ext = "pdf";return true;}
    if(a_format=="gl2ps_svg") {a_ext = "svg";return true;}
    if(a_format=="gl2ps_tex") {a_ext = "tex";return true;}
    if(a_format=="gl2ps_pgf") {a_ext = "pgf";return true;}
    if(a_format=="zb_ps")     {a_ext = "ps";return true;}
    if(a_format=="zb_png")    {a_ext = "png";return true;}
    if(a_format=="zb_jpeg")   {a_ext = "jpeg";return true;}
    a_ext.clear();
    return false;
  }

protected:
  std::string fFileName;
  std::string fFilePrefix;
  unsigned int fFileIndex;
  bool fResetFileIndex;
};

#endif
