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
//
//
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the OpenGL work by Andrew Walkden (27th March 1996).
// Vulkan viewer - opens window, hard copy, etc.

#ifndef G4VULKANVIEWER_HH
#define G4VULKANVIEWER_HH

#include "G4Vulkan.hh"
#include "G4VViewer.hh"

class G4VulkanSceneHandler;
class G4vk2ps;
class G4Text;

class G4VulkanViewerPickMap
{
  public:

    inline void setName(G4String n) { fName = n; }

    inline void setHitNumber(G4int n) { fHitNumber = n; }

    inline void setSubHitNumber(G4int n) { fSubHitNumber = n; }
    inline void setPickName(G4int n) { fPickName = n; }

    inline void addAttributes(G4String att) { fAttributes.push_back(att); }

    inline G4String getName() { return fName; }
    inline G4int getHitNumber() { return fHitNumber; }

    inline G4int getSubHitNumber() { return fSubHitNumber; }

    inline G4int getPickName() { return fPickName; }

    inline const std::vector<G4String>& getAttributes() const { return fAttributes; }

    G4String print();

  private:

    G4String fName;
    G4int fHitNumber;
    G4int fSubHitNumber;
    G4int fPickName;
    std::vector<G4String> fAttributes;
};

// Base class for various VulkanView classes.
class G4VulkanViewer : virtual public G4VViewer
{
    friend class G4VulkanSceneHandler;
    friend class G4VulkanImmediateSceneHandler;
    friend class G4VulkanStoredSceneHandler;
    friend class G4VulkanFileSceneHandler;
    friend class G4VulkanViewerMessenger;

  public:

    virtual void ClearView();
    void ClearViewWithoutFlush();

    virtual bool exportImage(std::string name = "", int width = -1, int height = -1);
    bool setExportImageFormat(std::string format, bool quiet = false);

  protected:

    G4VulkanViewer(G4VulkanSceneHandler& scene);
    virtual ~G4VulkanViewer();

  private:

    G4VulkanViewer(const G4VulkanViewer&);
    G4VulkanViewer& operator=(const G4VulkanViewer&);

  protected:

    virtual void SetView();
    virtual void ResetView();

    virtual void DrawText(const G4Text&);
    void ChangePointSize(G4double size);
    void ChangeLineWidth(G4double width);
    void HaloingFirstPass();
    void HaloingSecondPass();
    void HLRFirstPass();
    void HLRSecondPass();
    void HLRThirdPass();
    void InitializeVKView();
    void ResizeVKView();
    void ResizeWindow(unsigned int, unsigned int);
    virtual G4String Pick(VKdouble x, VKdouble y);
    const std::vector<G4VulkanViewerPickMap*>& GetPickDetails(VKdouble x, VKdouble y);
    virtual void CreateFontLists() {}
    void rotateScene(G4double dx, G4double dy);
    void rotateSceneToggle(G4double dx, G4double dy);
    //////////////////////////////Vectored PostScript production functions///
    // print EPS file. Depending of fVectoredPs, it will print Vectored or not
    void setExportSize(G4int, G4int);
    // set the new print size.
    // -1 means 'print size' = 'window size'
    // Setting size greater than max Vulkan viewport size will set the size to
    // maximum
    bool setExportFilename(G4String name, G4bool inc = true);
    // set export filename.
    // if inc, then the filename will be increment by one each time
    // try to guesss the correct format according to the extention

    std::string getRealPrintFilename();
    unsigned int getWinWidth() const;
    unsigned int getWinHeight() const;
    G4bool sizeHasChanged();
    // return true if size has change since last redraw
    G4double GetSceneNearWidth();
    VKdouble getSceneFarWidth();
    VKdouble getSceneDepth();
    void addExportImageFormat(std::string format);
    // add a image format to the available export format list
    G4bool isVk2psWriting();
    G4bool isFramebufferReady();

    void g4VkPickMatrix(VKdouble x, VKdouble y, VKdouble width, VKdouble height,
                         VKint viewport[4]);
    // MESA implementation of gluPickMatrix

    void g4VkLookAt(VKdouble eyex, VKdouble eyey, VKdouble eyez, VKdouble centerx,
                     VKdouble centery, VKdouble centerz, VKdouble upx, VKdouble upy, VKdouble upz);
    // MESA implementation of gluLookAt
    void g4VkOrtho(VKdouble left, VKdouble right, VKdouble bottom, VKdouble top, VKdouble near,
                   VKdouble far);
    // Redefinition on vkOrtho to solve precision issues
    void g4VkFrustum(VKdouble left, VKdouble right, VKdouble bottom, VKdouble top, VKdouble near,
                     VKdouble far);
    // Redefinition on vkFrustum to solve precision issues

    bool IsGettingPickInfos() const { return fIsGettingPickInfos; }

    G4bool fPrintColour;
    G4bool fVectoredPs;

    G4VulkanSceneHandler& fVulkanSceneHandler;
    G4Colour background;  // the Vulkan clear colour
    G4bool transparency_enabled,  // is alpha blending enabled?
      antialiasing_enabled,  // is antialiasing enabled?
      haloing_enabled;  // is haloing enabled for wireframe?
    G4vk2ps* fVK2PSAction;

    G4double fRot_sens;  // Rotation sensibility in degrees
    G4double fPan_sens;  // Translation sensibility
    unsigned int fWinSize_x;
    unsigned int fWinSize_y;
    std::vector<std::string> fExportImageFormatVector;
    std::string fDefaultExportImageFormat;
    std::string fExportImageFormat;
    int fExportFilenameIndex;
    G4int fPrintSizeX;
    G4int fPrintSizeY;

  private:

    G4float fPointSize;
    G4String fExportFilename;
    G4String fDefaultExportFilename;
    G4bool fSizeHasChanged;
    int fVk2psDefaultLineWith;
    int fVk2psDefaultPointSize;
    bool fVkViewInitialized;

    // size of the Vulkan frame
    void rotateSceneThetaPhi(G4double dx, G4double dy);
    void rotateSceneInViewDirection(G4double dx, G4double dy);
    bool printVk2PS();
    G4int getRealExportWidth();
    G4int getRealExportHeight();
    VKubyte* grabPixels(int inColor, unsigned int width, unsigned int height);
    bool printNonVectoredEPS();
    // print non vectored EPS files

    bool printVectoredEPS();
    // print vectored EPS files

    bool fIsGettingPickInfos;
    // Block SetView() during picking
};

#endif
