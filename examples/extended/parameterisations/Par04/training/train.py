"""
** train **
    - defines data loading parameters and calls the data preprocessing function 
    - defines the model parameters and instantiates the VAE model
    - performs the training
"""

# 1. Data loading/preprocessing
from utils import *
# Directory where the HDF5 files are saved
init_dir = './detector_'
# Number of calorimeter layers 
nCells_z = 45
# Segmentation in the r,phi direction
nCells_r = 18
nCells_phi = 50
# Total number of readout cells (represents the number of nodes in the input/output layers of the model)
original_dim = nCells_z*nCells_r*nCells_phi
# Minimum and maximum primary particle energy to consider for training in GeV units 
min_energy = 1
max_energy = 1024
# Minimum and maximum primary particle angle to consider for training in degrees units 
min_angle = 50
max_angle = 90
# The preprocess function reads the data and performs preprocessing and encoding for the values of energy, angle and geometry
energies_Train,condE_Train,condAngle_Train,condGeo_Train = preprocess(init_dir,original_dim,min_angle,max_angle,min_energy,max_energy)

# 2. Model architecture
import model
# Instantiate a VAE model and define all the parameters
vae = model.VAE(batch_size=100 ,
          original_dim=original_dim,
          intermediate_dim1=100,  
          intermediate_dim2=50,
          intermediate_dim3=20,
          intermediate_dim4=10+4,
          latent_dim=10,
          epsilon_std=1.,
          mu=0,
          epochs=10000,
          lr=0.001,
          activ=tf.keras.layers.LeakyReLU(),  
          outActiv='sigmoid',  
          validation_split=0.05,
          wReco=original_dim,
          wkl=0.5,
          optimizer=optimizers.Adam(),
          ki='RandomNormal',
          bi='Zeros',
          earlyStop=False,
          checkpoint_dir = "."
          )

# 3. Model training
history = vae.train(energies_Train,
        condE_Train,
        condAngle_Train,
        condGeo_Train
        )

# 4. Save the model (ony the decoder part) after traing 
self.vae.decoder.save("decoder.h5")

# 5. Convert the model to ONNX format
import keras2onnx
import tensorflow
# Create the Keras model and convert itinto an ONNX model
kerasModel = tensorflow.keras.models.load_model("decoder.h5") 
onnxModel = keras2onnx.convert_keras(kerasModel,"name")
# Save the ONNX model. Generator.onnx can then be used to perform the inference in the example
keras2onnx.save_model(onnxModel,"Generator.onnx")

"""
# In order to convert the model into a format that can be used with the LWTNN library
# 1. After training :
# serialize model to JSON
json_model = self.vae.decoder.to_json()
with open("decoder.json", "w") as json_file:
    json_file.write(json_model)
# serialize weights to HDF5
self.vae.decoder.save_weights("decoder.h5")
# 2. Externally, after building the LWTNN code available at https://github.com/lwtnn/lwtnn 
#   2.1 Run the kerasfunc2json python script (available in lwtnn/ converters/) to generate a template file of your functional model input variables by calling:
#       $ kerasfunc2json.py decoder.json decoder.h5 > inputs.json
#   2.2 Run again kerasfunc2json script to get your output file that would be used for the inference in the example
#       $ kerasfunc2json.py decoder.json decoder.h5 inputs.json > Generator.json
"""
