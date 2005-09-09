//
// File name:     RadmonTDetectorVolumesWithHoleDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithHoleDecorator.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with a hole a component with the following method
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH

 class RadmonTDetectorLayerConstructor;
 class RadmonDetectorLayerVolumesList;

 template <class LayerVolumesComponent> class RadmonTDetectorVolumesWithHoleDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithHoleDecorator(const RadmonTDetectorLayerConstructor * owner);
                                               ~RadmonTDetectorVolumesWithHoleDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithHoleDecorator();
                                                RadmonTDetectorVolumesWithHoleDecorator(const  RadmonTDetectorVolumesWithHoleDecorator & copy);
   RadmonTDetectorVolumesWithHoleDecorator &    operator=(const  RadmonTDetectorVolumesWithHoleDecorator & copy);

  // Private attributes
   const RadmonTDetectorLayerConstructor *      attributesOwner;
 };
#endif /* RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH */
