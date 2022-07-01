"""
** utils **
defines the data loading and preprocessing function 
"""

# Setup
import h5py
import numpy as np

# preprocess function returns the array of the shower energies and the condition arrays 
"""
    - init_dir: the name of the directory which contains the HDF5 files 
    - size_1DVec: represents the size of the input and output layer of the VAE which corresponds to the total number of readout cells
    - min_energy,max_energy: minimum and maximum primary particle energy to consider for training in GeV units  
    - min_angle and max_angle: minimum and maximum primary particle angle to consider for training in degrees units
"""
def preprocess(init_dir,size_1DVec,min_angle,max_angle,min_energy,max_energy):
    energies_Train = []
    condE_Train = []
    condAngle_Train = []
    condGeo_Train = []
    # This example is trained using 2 detector geometries
    for geo in [ 'SiW' , 'SciPb' ]:
        dirGeo = init_dir + geo + '/'
        energyParticle=min_energy
        # loop over the energies in powers of 2
        while(energyParticle<=max_energy):
            # loop over the angles in a step of 10
            for angleParticle in range(min_angle,max_angle+10,10):
                fName = 'Energy_%s_Angle_%s.hdf5' %(energyParticle,angleParticle)
                fName = dirGeo + fName
                # read the HDF5 file
                h5 = h5py.File(fName,'r')
                # get the key value of the group from the HDF5 file
                GroupKey = 'Grp_Angle_%s_E_%s'%(angleParticle,energyParticle)
                # get all key values of one group
                listKeys = list( h5[GroupKey].keys() )
                # loop over the events
                for ckey in listKeys:
                    # scale the energy of each cell to the energy of the primary particle (in MeV units) 
                    energyArray = np.array(h5[GroupKey][ckey])/(energyParticle*1000)
                    energies_Train.append( energyArray.reshape(size_1DVec)  )
                # build the energy and angle condition vectors
                condE_Train.append( [energyParticle/mamax_energyxE]*len(listKeys) )
                condAngle_Train.append( [angleParticle/max_angle]*len(listKeys) )
                # build the geometry condition vector (1 hot encoding vector)
                if( geo == 'SiW' ):
                    condGeo_Train.append( [[0,1]]*len(listKeys) )
                else:
                    condGeo_Train.append( [[1,0]]*len(listKeys) )
            energyParticle*=2
    # return numpy arrays 
    energies_Train = np.array(energies_Train)
    condE_Train = np.concatenate(condE_Train)
    condAngle_Train = np.concatenate(condAngle_Train)
    condGeo_Train = np.concatenate(condGeo_Train)
    return energies_Train,condE_Train,condAngle_Train,condGeo_Train 
