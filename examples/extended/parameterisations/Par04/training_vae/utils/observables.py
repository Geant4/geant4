from dataclasses import dataclass
from enum import Enum

import numpy as np

from core.constants import N_CELLS_Z, N_CELLS_R, SIZE_Z, SIZE_R


@dataclass
class Observable:
    """ An abstract class defining interface of all observables.

    Do not use this class directly.

    Attributes:
          _input: A numpy array with shape = (NE, R, PHI, Z), where NE stays for number of events.
    """
    _input: np.ndarray


class ProfileType(Enum):
    """ Enum class of various profile types.

    """
    LONGITUDINAL = 0
    LATERAL = 1


@dataclass
class Profile(Observable):
    """ An abstract class describing behaviour of LongitudinalProfile and LateralProfile.

    Do not use this class directly. Use LongitudinalProfile or LateralProfile instead.

    """

    def calc_profile(self) -> np.ndarray:
        pass

    def calc_first_moment(self) -> np.ndarray:
        pass

    def calc_second_moment(self) -> np.ndarray:
        pass


@dataclass
class LongitudinalProfile(Profile):
    """ A class defining observables related to LongitudinalProfile.

    Attributes:
        _energies_per_event: A numpy array with shape = (NE, Z) where NE stays for a number of events. An
            element [i, j] is a sum of energies detected in all cells located in a jth layer for an ith event.
        _total_energy_per_event: A numpy array with shape = (NE, ). An element [i] is a sum of energies detected in all
            cells for an ith event.
        _w: A numpy array = [0, 1, ..., Z - 1] which represents weights used in computation of first and second moment.

    """

    def __post_init__(self):
        self._energies_per_event = np.sum(self._input, axis=(1, 2))
        self._total_energy_per_event = np.sum(self._energies_per_event, axis=1)
        self._w = np.arange(N_CELLS_Z)

    def calc_profile(self) -> np.ndarray:
        """ Calculates a longitudinal profile.

        A longitudinal profile for a given layer l (l = 0, ..., Z - 1) is defined as:
        sum_{i = 0}^{NE - 1} energy_per_event[i, l].

        Returns:
            A numpy array of longitudinal profiles for each layer with a shape = (Z, ).

        """
        return np.sum(self._energies_per_event, axis=0)

    def calc_first_moment(self) -> np.ndarray:
        """ Calculates a first moment of profile.

        A first moment of a longitudinal profile for a given event e (e = 0, ..., NE - 1) is defined as:
        FM[e] = alpha * (sum_{i = 0}^{Z - 1} energies_per_event[e, i] * w[i]) / total_energy_per_event[e], where
        w = [0, 1, 2, ..., Z - 1],
        alpha = SIZE_Z defined in core/constants.py.

        Returns:
            A numpy array of first moments of longitudinal profiles for each event with a shape = (NE, ).

        """
        return SIZE_Z * np.dot(self._energies_per_event, self._w) / self._total_energy_per_event

    def calc_second_moment(self) -> np.ndarray:
        """ Calculates a second moment of a longitudinal profile.

        A second moment of a longitudinal profile for a given event e (e = 0, ..., NE - 1) is defined as:
        SM[e] = (sum_{i = 0}^{Z - 1} (w[i] - alpha - FM[e])^2 * energies_per_event[e, i]) total_energy_per_event[e],
        where
        w = [0, 1, 2, ..., Z - 1],
        alpha = SIZE_Z defined in ochre/constants.py

        Returns:
            A numpy array of second moments of longitudinal profiles for each event with a shape = (NE, ).
        """
        first_moment = self.calc_first_moment()
        first_moment = np.expand_dims(first_moment, axis=1)
        w = np.expand_dims(self._w, axis=0)
        # w has now a shape = [1, Z] and first moment has a shape = [NE, 1]. There is a broadcasting in the line
        # below how that one create an array with a shape = [NE, Z]
        return np.sum(np.multiply(np.power(w * SIZE_Z - first_moment, 2), self._energies_per_event),
                      axis=1) / self._total_energy_per_event


@dataclass
class LateralProfile(Profile):
    """ A class defining observables related to LateralProfile.

    Attributes:
        _energies_per_event: A numpy array with shape = (NE, R) where NE stays for a number of events. An
            element [i, j] is a sum of energies detected in all cells located in a jth layer for an ith event.
        _total_energy_per_event: A numpy array with shape = (NE, ). An element [i] is a sum of energies detected in all
            cells for an ith event.
        _w: A numpy array = [0, 1, ..., R - 1] which represents weights used in computation of first and second moment.

    """

    def __post_init__(self):
        self._energies_per_event = np.sum(self._input, axis=(2, 3))
        self._total_energy_per_event = np.sum(self._energies_per_event, axis=1)
        self._w = np.arange(N_CELLS_R)

    def calc_profile(self) -> np.ndarray:
        """ Calculates a lateral profile.

        A lateral profile for a given layer l (l = 0, ..., R - 1) is defined as:
        sum_{i = 0}^{NE - 1} energy_per_event[i, l].

        Returns:
            A numpy array of longitudinal profiles for each layer with a shape = (R, ).

        """
        return np.sum(self._energies_per_event, axis=0)

    def calc_first_moment(self) -> np.ndarray:
        """ Calculates a first moment of profile.

        A first moment of a lateral profile for a given event e (e = 0, ..., NE - 1) is defined as:
        FM[e] = alpha * (sum_{i = 0}^{R - 1} energies_per_event[e, i] * w[i]) / total_energy_per_event[e], where
        w = [0, 1, 2, ..., R - 1],
        alpha = SIZE_R defined in core/constants.py.

        Returns:
            A numpy array of first moments of lateral profiles for each event with a shape = (NE, ).

        """
        return SIZE_R * np.dot(self._energies_per_event, self._w) / self._total_energy_per_event

    def calc_second_moment(self) -> np.ndarray:
        """ Calculates a second moment of a lateral profile.

        A second moment of a lateral profile for a given event e (e = 0, ..., NE - 1) is defined as:
        SM[e] = (sum_{i = 0}^{R - 1} (w[i] - alpha - FM[e])^2 * energies_per_event[e, i]) total_energy_per_event[e],
        where
        w = [0, 1, 2, ..., R - 1],
        alpha = SIZE_R defined in ochre/constants.py

        Returns:
            A numpy array of second moments of lateral profiles for each event with a shape = (NE, ).
        """
        first_moment = self.calc_first_moment()
        first_moment = np.expand_dims(first_moment, axis=1)
        w = np.expand_dims(self._w, axis=0)
        # w has now a shape = [1, R] and first moment has a shape = [NE, 1]. There is a broadcasting in the line
        # below how that one create an array with a shape = [NE, R]
        return np.sum(np.multiply(np.power(w * SIZE_R - first_moment, 2), self._energies_per_event),
                      axis=1) / self._total_energy_per_event


@dataclass
class Energy(Observable):
    """ A class defining observables total energy per event and cell energy.

    """

    def calc_total_energy(self):
        """ Calculates total energy detected in an event.

        Total energy for a given event e (e = 0, ..., NE - 1) is defined as a sum of energies detected in all cells
        for this event.

        Returns:
            A numpy array of total energy values with shape = (NE, ).
        """
        return np.sum(self._input, axis=(1, 2, 3))

    def calc_cell_energy(self):
        """ Calculates cell energy.

        Cell energy for a given event (e = 0, ..., NE - 1) is defined by an array with shape (R * PHI * Z) storing
        values of energy in particular cells.

        Returns:
            A numpy array of cell energy values with shape = (NE * R * PHI * Z, ).

        """
        return np.copy(self._input).reshape(-1)

    def calc_energy_per_layer(self):
        """ Calculates total energy detected in a particular layer.

        Energy per layer for a given event (e = 0, ..., NE - 1) is defined by an array with shape (Z, ) storing
        values of total energy detected in a particular layer

        Returns:
            A numpy array of cell energy values with shape = (NE, Z).

        """
        return np.sum(self._input, axis=(1, 2))
