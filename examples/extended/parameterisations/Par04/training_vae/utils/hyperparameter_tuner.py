from dataclasses import dataclass
from typing import Tuple, Dict, Any, List

import numpy as np
from optuna import Trial, create_study, get_all_study_summaries, load_study
from optuna.pruners import MedianPruner
from optuna.samplers import TPESampler
from optuna.trial import TrialState

from core.constants import LEARNING_RATE, BATCH_SIZE_PER_REPLICA, ACTIVATION, OUT_ACTIVATION, \
    OPTIMIZER_TYPE, KERNEL_INITIALIZER, BIAS_INITIALIZER, N_TRIALS, LATENT_DIM, \
    INTERMEDIATE_DIMS, MAX_HIDDEN_LAYER_DIM, GLOBAL_CHECKPOINT_DIR
from core.model import VAEHandler
from utils.preprocess import preprocess


@dataclass
class HyperparameterTuner:
    """Tuner which looks for the best hyperparameters of a Variational Autoencoder specified in model.py.

    Currently, supported hyperparameters are: dimension of latent space, number of hidden layers, learning rate,
    activation function, activation function after the final layer, optimizer type, kernel initializer,
    bias initializer, batch size.

    Attributes:
        _discrete_parameters: A dictionary of hyperparameters taking discrete values in the range [low, high].
        _continuous_parameters: A dictionary of hyperparameters taking continuous values in the range [low, high].
        _categorical_parameters: A dictionary of hyperparameters taking values specified by the list of them.
        _storage: A string representing URL to a database required for a distributed training
        _study_name: A string, a name of study.

    """
    _discrete_parameters: Dict[str, Tuple[int, int]]
    _continuous_parameters: Dict[str, Tuple[float, float]]
    _categorical_parameters: Dict[str, List[Any]]
    _storage: str = None
    _study_name: str = None

    def _check_hyperparameters(self):
        available_hyperparameters = ["latent_dim", "nb_hidden_layers", "learning_rate", "activation", "out_activation",
                                     "optimizer_type", "kernel_initializer", "bias_initializer",
                                     "batch_size_per_replica"]
        hyperparameters_to_be_optimized = list(self._discrete_parameters.keys()) + list(
            self._continuous_parameters.keys()) + list(self._categorical_parameters.keys())
        for hyperparameter_name in hyperparameters_to_be_optimized:
            if hyperparameter_name not in available_hyperparameters:
                raise Exception(f"Unknown hyperparameter: {hyperparameter_name}")

    def __post_init__(self):
        self._check_hyperparameters()
        self._energies_train, self._cond_e_train, self._cond_angle_train, self._cond_geo_train = preprocess()

        if self._storage is not None and self._study_name is not None:
            # Parallel optimization
            study_summaries = get_all_study_summaries(self._storage)
            if any(self._study_name == study_summary.study_name for study_summary in study_summaries):
                # The study is already created in the database. Load it.
                self._study = load_study(self._study_name, self._storage)
            else:
                # The study does not exist in the database. Create a new one.
                self._study = create_study(storage=self._storage, sampler=TPESampler(), pruner=MedianPruner(),
                                           study_name=self._study_name, direction="minimize")
        else:
            # Single optimization
            self._study = create_study(sampler=TPESampler(), pruner=MedianPruner(), direction="minimize")

    def _create_model_handler(self, trial: Trial) -> VAEHandler:
        """For a given trail builds the model.

        Optuna suggests parameters like dimensions of particular layers of the model, learning rate, optimizer, etc.

        Args:
            trial: Optuna's trial

        Returns:
            Variational Autoencoder (VAE)
        """

        # Discrete parameters
        if "latent_dim" in self._discrete_parameters.keys():
            latent_dim = trial.suggest_int(name="latent_dim",
                                           low=self._discrete_parameters["latent_dim"][0],
                                           high=self._discrete_parameters["latent_dim"][1])
        else:
            latent_dim = LATENT_DIM

        if "nb_hidden_layers" in self._discrete_parameters.keys():
            nb_hidden_layers = trial.suggest_int(name="nb_hidden_layers",
                                                 low=self._discrete_parameters["nb_hidden_layers"][0],
                                                 high=self._discrete_parameters["nb_hidden_layers"][1])

            all_possible = np.arange(start=latent_dim + 5, stop=MAX_HIDDEN_LAYER_DIM)
            chunks = np.array_split(all_possible, nb_hidden_layers)
            ranges = [(chunk[0], chunk[-1]) for chunk in chunks]
            ranges = reversed(ranges)

            # Cast from np.int to int allows to become JSON serializable.
            intermediate_dims = [trial.suggest_int(name=f"intermediate_dim_{i}", low=int(low), high=int(high)) for
                                 i, (low, high)
                                 in enumerate(ranges)]
        else:
            intermediate_dims = INTERMEDIATE_DIMS

        if "batch_size_per_replica" in self._discrete_parameters.keys():
            batch_size_per_replica = trial.suggest_int(name="batch_size_per_replica",
                                                       low=self._discrete_parameters["batch_size_per_replica"][0],
                                                       high=self._discrete_parameters["batch_size_per_replica"][1])
        else:
            batch_size_per_replica = BATCH_SIZE_PER_REPLICA

        # Continuous parameters
        if "learning_rate" in self._continuous_parameters.keys():
            learning_rate = trial.suggest_float(name="learning_rate",
                                                low=self._continuous_parameters["learning_rate"][0],
                                                high=self._continuous_parameters["learning_rate"][1])
        else:
            learning_rate = LEARNING_RATE

        # Categorical parameters
        if "activation" in self._categorical_parameters.keys():
            activation = trial.suggest_categorical(name="activation",
                                                   choices=self._categorical_parameters["activation"])
        else:
            activation = ACTIVATION

        if "out_activation" in self._categorical_parameters.keys():
            out_activation = trial.suggest_categorical(name="out_activation",
                                                       choices=self._categorical_parameters["out_activation"])
        else:
            out_activation = OUT_ACTIVATION

        if "optimizer_type" in self._categorical_parameters.keys():
            optimizer_type = trial.suggest_categorical(name="optimizer_type",
                                                       choices=self._categorical_parameters["optimizer_type"])
        else:
            optimizer_type = OPTIMIZER_TYPE

        if "kernel_initializer" in self._categorical_parameters.keys():
            kernel_initializer = trial.suggest_categorical(name="kernel_initializer",
                                                           choices=self._categorical_parameters["kernel_initializer"])
        else:
            kernel_initializer = KERNEL_INITIALIZER

        if "bias_initializer" in self._categorical_parameters.keys():
            bias_initializer = trial.suggest_categorical(name="bias_initializer",
                                                         choices=self._categorical_parameters["bias_initializer"])
        else:
            bias_initializer = BIAS_INITIALIZER

        checkpoint_dir = f"{GLOBAL_CHECKPOINT_DIR}/{self._study_name}/trial_{trial.number:03d}"

        return VAEHandler(_wandb_project_name=self._study_name,
                          _wandb_tags=["hyperparameter tuning", f"trial {trial.number}"],
                          _batch_size_per_replica=batch_size_per_replica,
                          _intermediate_dims=intermediate_dims,
                          latent_dim=latent_dim,
                          _learning_rate=learning_rate,
                          _activation=activation,
                          _out_activation=out_activation,
                          _optimizer_type=optimizer_type,
                          _kernel_initializer=kernel_initializer,
                          _bias_initializer=bias_initializer,
                          _checkpoint_dir=checkpoint_dir,
                          _early_stop=True,
                          _save_model_every_epoch=False,
                          _save_best_model=True,
                          )

    def _objective(self, trial: Trial) -> float:
        """For a given trial trains the model and returns an average validation loss.

        Args:
            trial: Optuna's trial

        Returns: One float numer which is a validation loss. It can be either calculated as an average of k trainings
        performed in cross validation mode or is one number obtained from  validation on unseen before, some fraction
        of the dataset.
        """

        # Generate the trial model.
        model_handler = self._create_model_handler(trial)

        # Train the model.
        verbose = True
        histories = model_handler.train(self._energies_train, self._cond_e_train, self._cond_angle_train,
                                        self._cond_geo_train, verbose)

        # Return validation loss (currently it is treated as an objective goal). Notice that we take into account the
        # best model according to the validation loss.
        final_validation_losses = [np.min(history.history["val_loss"]) for history in histories]
        avg_validation_loss = np.mean(final_validation_losses).item()
        return avg_validation_loss

    def tune(self) -> None:
        """Main tuning function.

        Based on a given study, tunes the model and prints detailed information about the best trial (value of the
        objective function and adjusted parameters).
        """

        self._study.optimize(func=self._objective, n_trials=N_TRIALS, gc_after_trial=True)
        pruned_trials = self._study.get_trials(deepcopy=False, states=(TrialState.PRUNED,))
        complete_trials = self._study.get_trials(deepcopy=False, states=(TrialState.COMPLETE,))
        print("Study statistics: ")
        print("  Number of finished trials: ", len(self._study.trials))
        print("  Number of pruned trials: ", len(pruned_trials))
        print("  Number of complete trials: ", len(complete_trials))

        print("Best trial:")
        trial = self._study.best_trial

        print("  Value: ", trial.value)

        print("  Params: ")
        for key, value in trial.params.items():
            print(f"    {key}: {value}")
