"""
** convert **
defines the conversion function to and ONNX file
"""

import argparse
import sys

import tf2onnx
import numpy as np
from onnxruntime import InferenceSession

from core.constants import GLOBAL_CHECKPOINT_DIR, CONV_DIR, ORIGINAL_DIM
from core.model import VAEHandler
"""
    epoch: epoch of the saved checkpoint model
    study-name: study-name for which the model is trained for
"""


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--epoch", type=int, default=None)
    p.add_argument("--study-name", type=str, default="default_study_name")
    args = p.parse_args()
    return args


# main function
def main(argv):
    # 1. Set up the model to convert
    # Parse commandline arguments
    args = parse_args(argv)
    epoch = args.epoch
    study_name = args.study_name

    # Instantiate and load a saved model
    vae = VAEHandler()

    # Load the saved weights
    weights_dir = f"VAE_epoch_{epoch:03}" if epoch is not None else "VAE_best"
    vae.model.load_weights(
        f"{GLOBAL_CHECKPOINT_DIR}/{study_name}/{weights_dir}/model_weights"
    ).expect_partial()

    # 2. Convert the model to ONNX format
    # Create the Keras model, convert it into an ONNX model, and save.
    keras_model = vae.model.decoder
    output_path = f"{CONV_DIR}/{study_name}/Generator_{weights_dir}.onnx"
    onnx_model = tf2onnx.convert.from_keras(keras_model,
                                            output_path=output_path)

    # Checking the converted model
    input_1 = np.random.randn(10).astype(np.float32).reshape(1, -1)
    input_2 = np.random.randn(1).astype(np.float32).reshape(1, -1)
    input_3 = np.random.randn(1).astype(np.float32).reshape(1, -1)
    input_4 = np.random.randn(2).astype(np.float32).reshape(1, -1)

    sess = InferenceSession(output_path)
    # TODO: @Piyush-555 Find a way to use predefined names
    result = sess.run(
        None, {
            'input_9': input_1,
            'input_6': input_2,
            'input_7': input_3,
            'input_8': input_4
        })
    assert result[0].shape[1] == ORIGINAL_DIM


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
