from argparse import ArgumentParser

from core.constants import MAX_GPU_MEMORY_ALLOCATION, GPU_IDS
from utils.gpu_limiter import GPULimiter
from utils.optimizer import OptimizerType

# Hyperparemeters to be optimized.
discrete_parameters = {"nb_hidden_layers": (1, 6), "latent_dim": (15, 100)}
continuous_parameters = {"learning_rate": (0.0001, 0.005)}
categorical_parameters = {"optimizer_type": [OptimizerType.ADAM, OptimizerType.RMSPROP]}


def parse_args():
    argument_parser = ArgumentParser()
    argument_parser.add_argument("--study-name", type=str, default="default_study_name")
    argument_parser.add_argument("--storage", type=str)
    argument_parser.add_argument("--max-gpu-memory-allocation", type=int, default=MAX_GPU_MEMORY_ALLOCATION)
    argument_parser.add_argument("--gpu-ids", type=str, default=GPU_IDS)
    args = argument_parser.parse_args()
    return args


def main():
    # 0. Parse arguments.
    args = parse_args()
    study_name = args.study_name
    storage = args.storage
    max_gpu_memory_allocation = args.max_gpu_memory_allocation
    gpu_ids = args.gpu_ids

    # 1. Set GPU memory limits.
    GPULimiter(_gpu_ids=gpu_ids, _max_gpu_memory_allocation=max_gpu_memory_allocation)()

    # 2. Manufacture hyperparameter tuner.

    # This import must be local because otherwise it is impossible to call GPULimiter.
    from utils.hyperparameter_tuner import HyperparameterTuner
    hyperparameter_tuner = HyperparameterTuner(discrete_parameters, continuous_parameters, categorical_parameters,
                                               storage, study_name)

    # 3. Run main tuning function.
    hyperparameter_tuner.tune()
    # Watch out! This script neither deletes the study in DB nor deletes the database itself. If you are using
    # parallelized optimization, then you should care about deleting study in the database by yourself.


if __name__ == "__main__":
    exit(main())
