## CaloDiT-2

CaloDiT-2 is transformer-based diffusion model, which can be easily adapted to a new detector geometry.
This Par04 repository contains the ONNX and TorchScript versions of the model and how to use them for
Par04-SiW detector. *Note that the virtual cylindrical mesh is less granular than VAE.*

The source code for CaloDiT-2 can be found [here](https://gitlab.cern.ch/fastsim/diffusion4sim/-/tree/CaloDiT_v1?ref_type=tags).
It contains the training, adaptation, and distillation scripts along with the pretrained models. You can download the
pretrained models (not the .onnx/.pt files available with this repository), which acts as a checkpoint and finetune
it on the new dataset.

You can also modify the virtual mesh size and hence the architecture if needed. But in that case, you won't be able
to use the pretrained models. Pretrained models adopt a mesh as in [CaloChallenge Dataset-2](https://calochallenge.github.io/homepage/).
