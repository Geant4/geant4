//
// File name:     RadmonTDetectorCarvedVolumesDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorCarvedVolumesDecorator.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Carves the borders of a component with the following method
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH
 #define  RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH

 // Forward declaration
 class RadmonTDetectorLayerConstructor;
 class RadmonDetectorLayerVolumesList;

 template <class LayerVolumesComponent> class RadmonTDetectorCarvedVolumesDecorator
 {
  public:
                                                RadmonTDetectorCarvedVolumesDecorator(const RadmonTDetectorLayerConstructor * owner);
                                               ~RadmonTDetectorCarvedVolumesDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorCarvedVolumesDecorator();
                                                RadmonTDetectorCarvedVolumesDecorator(const RadmonTDetectorCarvedVolumesDecorator & copy);
   RadmonTDetectorCarvedVolumesDecorator &      operator=(const RadmonTDetectorCarvedVolumesDecorator & copy);

  // Private attributes
   const RadmonTDetectorLayerConstructor *      attributesOwner;
 };
#endif /* RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH */
