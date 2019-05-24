===================
PowerGrid Utilities
===================

Utilities for working with PowerGrid


.. code-block:: MATLAB

   % Load files from PG
   cd TO_DIRECTORY_WITH_PG_OUTPUT
   img = collectPowerGridImgOutput();

   % Reshuffle multiband slabs into 3D volumes
   img3D = reshuffleMultibandPowerGridOutput(img);
