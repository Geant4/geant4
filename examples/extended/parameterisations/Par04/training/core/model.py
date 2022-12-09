import gc
from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np
import tensorflow as tf
import wandb
from sklearn.model_selection import KFold
from tensorflow.keras import backend as K
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, History, Callback
from tensorflow.keras.layers import BatchNormalization, Input, Dense, Layer, concatenate
from tensorflow.keras.losses import BinaryCrossentropy, Reduction
from tensorflow.keras.models import Model
from tensorflow.python.data import Dataset
from tensorflow.python.distribute.distribute_lib import Strategy
from tensorflow.python.distribute.mirrored_strategy import MirroredStrategy
from wandb.keras import WandbCallback

from core.constants import ORIGINAL_DIM, LATENT_DIM, BATCH_SIZE_PER_REPLICA, EPOCHS, LEARNING_RATE, ACTIVATION, \
    OUT_ACTIVATION, OPTIMIZER_TYPE, KERNEL_INITIALIZER, GLOBAL_CHECKPOINT_DIR, EARLY_STOP, BIAS_INITIALIZER, \
    INTERMEDIATE_DIMS, SAVE_MODEL_EVERY_EPOCH, SAVE_BEST_MODEL, PATIENCE, MIN_DELTA, BEST_MODEL_FILENAME, \
    NUMBER_OF_K_FOLD_SPLITS, VALIDATION_SPLIT, WANDB_ENTITY
from utils.optimizer import OptimizerFactory, OptimizerType


class _Sampling(Layer):
    """ Custom layer to do the reparameterization trick: sample random latent vectors z from the latent Gaussian
    distribution.

    The sampled vector z is given by sampled_z = mean + std * epsilon
    """

    def __call__(self, inputs, **kwargs):
        z_mean, z_log_var, epsilon = inputs
        z_sigma = K.exp(0.5 * z_log_var)
        return z_mean + z_sigma * epsilon


# KL divergence computation
class _KLDivergenceLayer(Layer):

    def call(self, inputs, **kwargs):
        mu, log_var = inputs
        kl_loss = -0.5 * (1 + log_var - K.square(mu) - K.exp(log_var))
        kl_loss = K.mean(K.sum(kl_loss, axis=-1))
        self.add_loss(kl_loss)
        return inputs


class VAE(Model):
    def get_config(self):
        config = super().get_config()
        config["encoder"] = self.encoder
        config["decoder"] = self.decoder
        return config

    def call(self, inputs, training=None, mask=None):
        _, e_input, angle_input, geo_input, _ = inputs
        z = self.encoder(inputs)
        return self.decoder([z, e_input, angle_input, geo_input])

    def __init__(self, encoder, decoder, **kwargs):
        super(VAE, self).__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        self._set_inputs(inputs=self.encoder.inputs, outputs=self(self.encoder.inputs))


@dataclass
class VAEHandler:
    """
    Class to handle building and training VAE models.
    """
    _wandb_project_name: str = None
    _wandb_tags: List[str] = field(default_factory=list)
    _original_dim: int = ORIGINAL_DIM
    latent_dim: int = LATENT_DIM
    _batch_size_per_replica: int = BATCH_SIZE_PER_REPLICA
    _intermediate_dims: List[int] = field(default_factory=lambda: INTERMEDIATE_DIMS)
    _learning_rate: float = LEARNING_RATE
    _epochs: int = EPOCHS
    _activation: str = ACTIVATION
    _out_activation: str = OUT_ACTIVATION
    _number_of_k_fold_splits: float = NUMBER_OF_K_FOLD_SPLITS
    _optimizer_type: OptimizerType = OPTIMIZER_TYPE
    _kernel_initializer: str = KERNEL_INITIALIZER
    _bias_initializer: str = BIAS_INITIALIZER
    _checkpoint_dir: str = GLOBAL_CHECKPOINT_DIR
    _early_stop: bool = EARLY_STOP
    _save_model_every_epoch: bool = SAVE_MODEL_EVERY_EPOCH
    _save_best_model: bool = SAVE_BEST_MODEL
    _patience: int = PATIENCE
    _min_delta: float = MIN_DELTA
    _best_model_filename: str = BEST_MODEL_FILENAME
    _validation_split: float = VALIDATION_SPLIT
    _strategy: Strategy = MirroredStrategy()

    def __post_init__(self) -> None:
        # Calculate true batch size.
        self._batch_size = self._batch_size_per_replica * self._strategy.num_replicas_in_sync
        self._build_and_compile_new_model()
        # Setup Wandb.
        if self._wandb_project_name is not None:
            self._setup_wandb()

    def _setup_wandb(self) -> None:
        config = {
            "learning_rate": self._learning_rate,
            "batch_size": self._batch_size,
            "epochs": self._epochs,
            "optimizer_type": self._optimizer_type,
            "intermediate_dims": self._intermediate_dims,
            "latent_dim": self.latent_dim
        }
        # Reinit flag is needed for hyperparameter tuning. Whenever new training is started, new Wandb run should be
        # created.
        wandb.init(project=self._wandb_project_name, entity=WANDB_ENTITY, reinit=True, config=config,
                   tags=self._wandb_tags)

    def _build_and_compile_new_model(self) -> None:
        """ Builds and compiles a new model.

        VAEHandler keep a list of VAE instance. The reason is that while k-fold cross validation is performed,
        each fold requires a new, clear instance of model. New model is always added at the end of the list of
        existing ones.

        Returns: None

        """
        # Build encoder and decoder.
        encoder = self._build_encoder()
        decoder = self._build_decoder()

        # Compile model within a distributed strategy.
        with self._strategy.scope():
            # Build VAE.
            self.model = VAE(encoder, decoder)
            # Manufacture an optimizer and compile model with.
            optimizer = OptimizerFactory.create_optimizer(self._optimizer_type, self._learning_rate)
            reconstruction_loss = BinaryCrossentropy(reduction=Reduction.SUM)
            self.model.compile(optimizer=optimizer, loss=[reconstruction_loss], loss_weights=[ORIGINAL_DIM])

    def _prepare_input_layers(self, for_encoder: bool) -> List[Input]:
        """
        Create four Input layers. Each of them is responsible to take respectively: batch of showers/batch of latent
        vectors, batch of energies, batch of angles, batch of geometries.

        Args:
            for_encoder: Boolean which decides whether an input is full dimensional shower or a latent vector.

        Returns:
            List of Input layers (five for encoder and four for decoder).

        """
        e_input = Input(shape=(1,))
        angle_input = Input(shape=(1,))
        geo_input = Input(shape=(2,))
        if for_encoder:
            x_input = Input(shape=self._original_dim)
            eps_input = Input(shape=self.latent_dim)
            return [x_input, e_input, angle_input, geo_input, eps_input]
        else:
            x_input = Input(shape=self.latent_dim)
            return [x_input, e_input, angle_input, geo_input]

    def _build_encoder(self) -> Model:
        """ Based on a list of intermediate dimensions, activation function and initializers for kernel and bias builds
        the encoder.

        Returns:
             Encoder is returned as a keras.Model.

        """

        with self._strategy.scope():
            # Prepare input layer.
            x_input, e_input, angle_input, geo_input, eps_input = self._prepare_input_layers(for_encoder=True)
            x = concatenate([x_input, e_input, angle_input, geo_input])
            # Construct hidden layers (Dense and Batch Normalization).
            for intermediate_dim in self._intermediate_dims:
                x = Dense(units=intermediate_dim, activation=self._activation,
                          kernel_initializer=self._kernel_initializer,
                          bias_initializer=self._bias_initializer)(x)
                x = BatchNormalization()(x)
            # Add Dense layer to get description of multidimensional Gaussian distribution in terms of mean
            # and log(variance).
            z_mean = Dense(self.latent_dim, name="z_mean")(x)
            z_log_var = Dense(self.latent_dim, name="z_log_var")(x)
            # Add KLDivergenceLayer responsible for calculation of KL loss.
            z_mean, z_log_var = _KLDivergenceLayer()([z_mean, z_log_var])
            # Sample a probe from the distribution.
            encoder_output = _Sampling()([z_mean, z_log_var, eps_input])
            # Create model.
            encoder = Model(inputs=[x_input, e_input, angle_input, geo_input, eps_input], outputs=encoder_output,
                            name="encoder")
        return encoder

    def _build_decoder(self) -> Model:
        """ Based on a list of intermediate dimensions, activation function and initializers for kernel and bias builds
        the decoder.

        Returns:
             Decoder is returned as a keras.Model.

        """

        with self._strategy.scope():
            # Prepare input layer.
            latent_input, e_input, angle_input, geo_input = self._prepare_input_layers(for_encoder=False)
            x = concatenate([latent_input, e_input, angle_input, geo_input])
            # Construct hidden layers (Dense and Batch Normalization).
            for intermediate_dim in reversed(self._intermediate_dims):
                x = Dense(units=intermediate_dim, activation=self._activation,
                          kernel_initializer=self._kernel_initializer,
                          bias_initializer=self._bias_initializer)(x)
                x = BatchNormalization()(x)
            # Add Dense layer to get output which shape is compatible in an input's shape.
            decoder_outputs = Dense(units=self._original_dim, activation=self._out_activation)(x)
            # Create model.
            decoder = Model(inputs=[latent_input, e_input, angle_input, geo_input], outputs=decoder_outputs,
                            name="decoder")
        return decoder

    def _manufacture_callbacks(self) -> List[Callback]:
        """
        Based on parameters set by the user, manufacture callbacks required for training.

        Returns:
            A list of `Callback` objects.

        """
        callbacks = []
        # If the early stopping flag is on then stop the training when a monitored metric (validation) has stopped
        # improving after (patience) number of epochs.
        if self._early_stop:
            callbacks.append(
                EarlyStopping(monitor="val_loss",
                              min_delta=self._min_delta,
                              patience=self._patience,
                              verbose=True,
                              restore_best_weights=True))
        # Save model after every epoch.
        if self._save_model_every_epoch:
            callbacks.append(ModelCheckpoint(filepath=f"{self._checkpoint_dir}/VAE_epoch_{{epoch:03}}/model_weights",
                                             monitor="val_loss",
                                             verbose=True,
                                             save_weights_only=True,
                                             mode="min",
                                             save_freq="epoch"))
        # Pass metadata to wandb.
        callbacks.append(WandbCallback(
            monitor="val_loss", verbose=0, mode="auto", save_model=False))
        return callbacks

    def _get_train_and_val_data(self, dataset: np.array, e_cond: np.array, angle_cond: np.array, geo_cond: np.array,
                                noise: np.array, train_indexes: np.array, validation_indexes: np.array) \
            -> Tuple[Dataset, Dataset]:
        """
        Splits data into train and validation set based on given lists of indexes.

        """

        # Prepare training data.
        train_dataset = dataset[train_indexes, :]
        train_e_cond = e_cond[train_indexes]
        train_angle_cond = angle_cond[train_indexes]
        train_geo_cond = geo_cond[train_indexes, :]
        train_noise = noise[train_indexes, :]

        # Prepare validation data.
        val_dataset = dataset[validation_indexes, :]
        val_e_cond = e_cond[validation_indexes]
        val_angle_cond = angle_cond[validation_indexes]
        val_geo_cond = geo_cond[validation_indexes, :]
        val_noise = noise[validation_indexes, :]

        # Gather them into tuples.
        train_x = (train_dataset, train_e_cond, train_angle_cond, train_geo_cond, train_noise)
        train_y = train_dataset
        val_x = (val_dataset, val_e_cond, val_angle_cond, val_geo_cond, val_noise)
        val_y = val_dataset

        # Wrap data in Dataset objects.
        # TODO(@mdragula): This approach requires loading the whole data set to RAM. It
        #  would be better to read the data partially when needed. Also one should bare in mind that using tf.Dataset
        #  slows down training process.
        train_data = Dataset.from_tensor_slices((train_x, train_y))
        val_data = Dataset.from_tensor_slices((val_x, val_y))

        # The batch size must now be set on the Dataset objects.
        train_data = train_data.batch(self._batch_size)
        val_data = val_data.batch(self._batch_size)

        # Disable AutoShard.
        options = tf.data.Options()
        options.experimental_distribute.auto_shard_policy = tf.data.experimental.AutoShardPolicy.DATA
        train_data = train_data.with_options(options)
        val_data = val_data.with_options(options)

        return train_data, val_data

    def _k_fold_training(self, dataset: np.array, e_cond: np.array, angle_cond: np.array, geo_cond: np.array,
                         noise: np.array, callbacks: List[Callback], verbose: bool = True) -> List[History]:
        """
        Performs K-fold cross validation training.

        Number of fold is defined by (self._number_of_k_fold_splits). Always shuffle the dataset.

        Args:
            dataset: A matrix representing showers. Shape =
                (number of samples, ORIGINAL_DIM = N_CELLS_Z * N_CELLS_R * N_CELLS_PHI).
            e_cond: A matrix representing an energy for each sample. Shape = (number of samples, ).
            angle_cond: A matrix representing an angle for each sample. Shape = (number of samples, ).
            geo_cond: A matrix representing a geometry of the detector for each sample. Shape = (number of samples, 2).
            noise: A matrix representing an additional noise needed to perform a reparametrization trick.
            callbacks: A list of callback forwarded to the fitting function.
            verbose: A boolean which says there the training should be performed in a verbose mode or not.

        Returns: A list of `History` objects.`History.history` attribute is a record of training loss values and
        metrics values at successive epochs, as well as validation loss values and validation metrics values (if
        applicable).

        """
        # TODO(@mdragula): KFold cross validation can be parallelized. Each fold is independent from each the others.
        k_fold = KFold(n_splits=self._number_of_k_fold_splits, shuffle=True)
        histories = []

        for i, (train_indexes, validation_indexes) in enumerate(k_fold.split(dataset)):
            print(f"K-fold: {i + 1}/{self._number_of_k_fold_splits}...")
            train_data, val_data = self._get_train_and_val_data(dataset, e_cond, angle_cond, geo_cond, noise,
                                                                train_indexes, validation_indexes)

            self._build_and_compile_new_model()

            history = self.model.fit(x=train_data,
                                     shuffle=True,
                                     epochs=self._epochs,
                                     verbose=verbose,
                                     validation_data=val_data,
                                     callbacks=callbacks
                                     )
            histories.append(history)

            if self._save_best_model:
                self.model.save_weights(f"{self._checkpoint_dir}/VAE_fold_{i + 1}/model_weights")
                print(f"Best model from fold {i + 1} was saved.")

            # Remove all unnecessary data from previous fold.
            del self.model
            del train_data
            del val_data
            tf.keras.backend.clear_session()
            gc.collect()

        return histories

    def _single_training(self, dataset: np.array, e_cond: np.array, angle_cond: np.array, geo_cond: np.array,
                         noise: np.ndarray, callbacks: List[Callback], verbose: bool = True) -> List[History]:
        """
        Performs a single training.

        A fraction of dataset (self._validation_split) is used as a validation data.

        Args:
            dataset: A matrix representing showers. Shape =
                (number of samples, ORIGINAL_DIM = N_CELLS_Z * N_CELLS_R * N_CELLS_PHI).
            e_cond: A matrix representing an energy for each sample. Shape = (number of samples, ).
            angle_cond: A matrix representing an angle for each sample. Shape = (number of samples, ).
            geo_cond: A matrix representing a geometry of the detector for each sample. Shape = (number of samples, 2).
            noise: A matrix representing an additional noise needed to perform a reparametrization trick.
            callbacks: A list of callback forwarded to the fitting function.
            verbose: A boolean which says there the training should be performed in a verbose mode or not.

        Returns: A one-element list of `History` objects.`History.history` attribute is a record of training loss
        values and metrics values at successive epochs, as well as validation loss values and validation metrics
        values (if applicable).

        """
        dataset_size, _ = dataset.shape
        permutation = np.random.permutation(dataset_size)
        split = int(dataset_size * self._validation_split)
        train_indexes, validation_indexes = permutation[split:], permutation[:split]

        train_data, val_data = self._get_train_and_val_data(dataset, e_cond, angle_cond, geo_cond, noise, train_indexes,
                                                            validation_indexes)

        history = self.model.fit(x=train_data,
                                 shuffle=True,
                                 epochs=self._epochs,
                                 verbose=verbose,
                                 validation_data=val_data,
                                 callbacks=callbacks
                                 )
        if self._save_best_model:
            self.model.save_weights(f"{self._checkpoint_dir}/VAE_best/model_weights")
            print("Best model was saved.")

        return [history]

    def train(self, dataset: np.array, e_cond: np.array, angle_cond: np.array, geo_cond: np.array,
              verbose: bool = True) -> List[History]:
        """
        For a given input data trains and validates the model.

        If the numer of K-fold splits > 1 then it runs K-fold cross validation, otherwise it runs a single training
        which uses (self._validation_split * 100) % of dataset as a validation data.

        Args:
            dataset: A matrix representing showers. Shape =
                (number of samples, ORIGINAL_DIM = N_CELLS_Z * N_CELLS_R * N_CELLS_PHI).
            e_cond: A matrix representing an energy for each sample. Shape = (number of samples, ).
            angle_cond: A matrix representing an angle for each sample. Shape = (number of samples, ).
            geo_cond: A matrix representing a geometry of the detector for each sample. Shape = (number of samples, 2).
            verbose: A boolean which says there the training should be performed in a verbose mode or not.

        Returns: A list of `History` objects.`History.history` attribute is a record of training loss values and
        metrics values at successive epochs, as well as validation loss values and validation metrics values (if
        applicable).

        """

        callbacks = self._manufacture_callbacks()

        noise = np.random.normal(0, 1, size=(dataset.shape[0], self.latent_dim))

        if self._number_of_k_fold_splits > 1:
            return self._k_fold_training(dataset, e_cond, angle_cond, geo_cond, noise, callbacks, verbose)
        else:
            return self._single_training(dataset, e_cond, angle_cond, geo_cond, noise, callbacks, verbose)
