"""
** setup **
creates necessary folders
"""

import os

from core.constants import INIT_DIR, GLOBAL_CHECKPOINT_DIR, CONV_DIR, VALID_DIR, GEN_DIR

for folder in [INIT_DIR,  # Directory to load the full simulation dataset
               GLOBAL_CHECKPOINT_DIR,  # Directory to save VAE checkpoints
               CONV_DIR,  # Directory to save model after conversion to a format that can be used in C++
               VALID_DIR,  # Directory to save validation plots
               GEN_DIR,  # Directory to save VAE generated showers
               ]:
    os.system(f"mkdir {folder}")
