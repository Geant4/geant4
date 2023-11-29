This repository contains the set of scripts used to train, generate and validate the generative model used
in this example.

- core/constants.py: defines the set of common variables.
- core/model.py: defines the VAE model class and a handler to construct the model.
- utils/preprocess.py: defines the data loading and preprocessing functions.
- utils/hyperparameter_tuner.py: defines the HyperparameterTuner class.
- utils/gpu_limiter.py: defines a logic responsible for GPU memory management.
- utils/observables.py: defines a set of observable possibly calculated from a shower.
- utils/plotter.py: defines plotting classes responsible for manufacturing various plots of observables.
- train.py: performs model training.
- generate.py: generate showers using a saved VAE model.
- observables.py: defines a set of shower observables.
- validate.py: creates validation plots using shower observables.
- convert.py: defines the conversion function to an ONNX file.
- tune_model.py: performs hyperparameters optimization.

## Getting Started

`setup.py` script creates necessary folders used to save model checkpoints, generate showers and validation plots.

```
python3 setup.py
``` 

## Full simulation dataset

The full simulation dataset can be downloaded from/linked to [Zenodo](https://zenodo.org/record/6082201#.Ypo5UeDRaL4).

## Training

In order to launch the training:

```
python3 train.py
``` 

You may specify those three following flags. If you do not, then default values will be used.

```--max-gpu-memory-allocation``` specifies a maximum memory allocation on a single, logic GPU unit. Should be given as
an integer.

```--gpu-ids``` specifies IDs of physical GPUs. Should be given as a string, separated with comas, no spaces.
If you specify more than one GPU then automatically ```tf.distribute.MirroredStrategy``` will be applied to the
training.

```--study-name``` specifies a study name. This name is used as an experiment name in W&B dashboard and as a name of
directory for saving models.

## Hyperparameters tuning

If you want to tune hyperparameters, specify in `tune_model.py` parameters to be tuned. There are three types of
parameters: discrete, continuous and categorical. Discrete and continuous require range specification (low, high), while
the categorical parameter requires a list of possible values to be chosen. Then run it with:

```
python3 tune_model.py
```

If you want to parallelize tuning process you need to specify a common storage (preferable MySQL database) by
setting `--storage="URL_TO_MYSQL_DATABASE"`. Then you can run multiple processes with the same command:

```
python3 tune_model.py --storage="URL_TO_MYSQL_DATABASE"
```

Similarly to training procedure, you may specify ```--max-gpu-memory-allocation```, ```--gpu-ids``` and
```--study-name```.

## ML shower generation (MLFastSim)

In order to generate showers using the ML model, use `generate.py` script and specify information of geometry, energy
and angle of the particle and the epoch of the saved checkpoint model. The number of events to generate can also be
specified (by default is set to 10.000):

```
python3 generate.py --geometry=SiW --energy=64 --angle=90 --epoch=1000 --study-name=YOUR_STUDY_NAME
``` 

If you do not specify an epoch number the based model (saved as ```VAEbest```) will be used for shower generation.

## Validation

In order to validate the MLFastSim and the full simulation, use `validate.py` script and specify information of
geometry, energy and angle of the particle:

```
python3 validate.py --geometry=SiW --energye=64 --angle=90 
``` 

## Conversion

After training and validation, the model can be converted into a format that can be used in C++, such as ONNX,
use `convert.py` script:

```
python3 convert.py --epoch 1000
```
