from enum import IntEnum

from tensorflow.keras.optimizers import Optimizer, Adadelta, Adagrad, Adam, Adamax, Ftrl, SGD, Nadam, RMSprop


class OptimizerType(IntEnum):
    """ Enum class of various optimizer types.

    This class must be IntEnum to be JSON serializable. This feature is important because, when Optuna's study is
    saved in a relational DB, all objects must be JSON serializable.
    """

    SGD = 0
    RMSPROP = 1
    ADAM = 2
    ADADELTA = 3
    ADAGRAD = 4
    ADAMAX = 5
    NADAM = 6
    FTRL = 7


class OptimizerFactory:
    """Factory of optimizer like Stochastic Gradient Descent, RMSProp, Adam, etc.
    """

    @staticmethod
    def create_optimizer(optimizer_type: OptimizerType, learning_rate: float) -> Optimizer:
        """For a given type and a learning rate creates an instance of optimizer.

        Args:
            optimizer_type: a type of optimizer
            learning_rate: a learning rate that should be passed to an optimizer

        Returns:
            An instance of optimizer.

        """
        if optimizer_type == OptimizerType.SGD:
            return SGD(learning_rate)
        elif optimizer_type == OptimizerType.RMSPROP:
            return RMSprop(learning_rate)
        elif optimizer_type == OptimizerType.ADAM:
            return Adam(learning_rate)
        elif optimizer_type == OptimizerType.ADADELTA:
            return Adadelta(learning_rate)
        elif optimizer_type == OptimizerType.ADAGRAD:
            return Adagrad(learning_rate)
        elif optimizer_type == OptimizerType.ADAMAX:
            return Adamax(learning_rate)
        elif optimizer_type == OptimizerType.NADAM:
            return Nadam(learning_rate)
        else:
            # i.e. optimizer_type == OptimizerType.FTRL
            return Ftrl(learning_rate)
