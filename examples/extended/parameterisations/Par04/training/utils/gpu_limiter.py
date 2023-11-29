import os
from dataclasses import dataclass

import tensorflow as tf


@dataclass
class GPULimiter:
    """
    Class responsible to set the limits of possible GPU usage by TensorFlow. Currently, the limiter creates one
    instance of logical device per physical device. This can be changed in a future.

    Attributes:
        _gpu_ids: A string representing visible devices for the process. Identifiers of physical GPUs should
            be separated by commas (no spaces).
        _max_gpu_memory_allocation: An integer specifying limit of allocated memory per logical device.

    """
    _gpu_ids: str
    _max_gpu_memory_allocation: int

    def __call__(self):
        os.environ["CUDA_VISIBLE_DEVICES"] = f"{self._gpu_ids}"
        gpus = tf.config.list_physical_devices('GPU')
        if gpus:
            # Restrict TensorFlow to only allocate max_gpu_memory_allocation*1024 MB of memory on one of the GPUs
            try:
                for gpu in gpus:
                    tf.config.set_logical_device_configuration(
                        gpu,
                        [tf.config.LogicalDeviceConfiguration(memory_limit=1024 * self._max_gpu_memory_allocation)])
                logical_gpus = tf.config.list_logical_devices('GPU')
                print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
            except RuntimeError as e:
                # Virtual devices must be set before GPUs have been initialized
                print(e)
