============
Installation
============


Git repository
--------------

Unfortunately, remo is not publicly available on github due to licensing issues of the
embedded legacy code. However, if you have access to the GERICS gitlab, you can get a copy
of the model code by simply cloning the repository:

.. code-block:: bash

   git clone https://git.gerics.de/REMO/remo2.git
   
to a directory of your choice.


Choosing a configuration
------------------------

Remo comes with a number of different model components and configurations that you can choose. This is usefull because
a regional model should be flexible enough to be easily adapted to the climatology of the region of interest. However,
to create a basic model setup, you can simply choose a standard configuration which is mostly similar to older remo version:

.. code-block:: bash

   ./setup -auto Remo2015 -objdir=Remo2015
