``io``
======
Interfacing of Tudat(py) to and from other applications.


This module provides with different functionalities that allow to export results from Tudat(py)
to other softwares such as MATLAB, to post-process results. A set of methods are also provided
to read inputs from other softwares and integrate them with Tudat(py).












Functions
---------
.. currentmodule:: tudatpy.io

.. autosummary::

   save2txt

   save_time_history_to_file

   get_resource_path

   get_ephemeris_path

   get_earth_orientation_path

   get_quadrature_path

   get_spice_kernel_path

   get_atmosphere_tables_path

   get_gravity_models_path

   get_space_weather_path

   read_vector_history_from_file

   read_matrix_history_from_file



.. autofunction:: tudatpy.io.save2txt

.. autofunction:: tudatpy.io.save_time_history_to_file

.. autofunction:: tudatpy.io.get_resource_path

.. autofunction:: tudatpy.io.get_ephemeris_path

.. autofunction:: tudatpy.io.get_earth_orientation_path

.. autofunction:: tudatpy.io.get_quadrature_path

.. autofunction:: tudatpy.io.get_spice_kernel_path

.. autofunction:: tudatpy.io.get_atmosphere_tables_path

.. autofunction:: tudatpy.io.get_gravity_models_path

.. autofunction:: tudatpy.io.get_space_weather_path

.. autofunction:: tudatpy.io.read_vector_history_from_file

.. autofunction:: tudatpy.io.read_matrix_history_from_file




Enumerations
------------
.. currentmodule:: tudatpy.io

.. autosummary::

   StaticCoefficientNames

   DynamicCoefficientNames



.. autoclass:: tudatpy.io.StaticCoefficientNames
   :members:

.. autoclass:: tudatpy.io.DynamicCoefficientNames
   :members:




Classes
-------
.. currentmodule:: tudatpy.io

.. autosummary::

   missile_DATCOM_data



.. autoclass:: tudatpy.io.missile_DATCOM_data
   :members:
   :special-members: __init__





