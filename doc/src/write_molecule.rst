.. index:: write_molecule

write_molecule command
======================

Syntax
""""""

.. code-block:: LAMMPS

   write_dump mol-ID file

* mol-ID = ID of the molecule template to be written
* file = name of file to write the molecule template to

Examples
""""""""

.. code-block:: LAMMPS

   write_molecule mol1 molecule1.mol
   write_molecule mol1 molecule1.json
   write_molecule twomols template_set%.mol

Description
"""""""""""

.. versionadded:: TBD

Write the data from a :doc:`molecule template<molecule>` to a file.

The file format is determined by the file name: if the file name ends
in ``.json`` the file will be written `JSON format <https://www.json.org/>`_,
otherwise the file is written in the native LAMMPS file format.

When the molecule template contains multiple data sets, the filename
must contain a '%' character which will be replaced by the data set number
so that each data set can be written to a separate file.

----------

Restrictions
""""""""""""

None

Related commands
""""""""""""""""

:doc:`molecule <molecule>`

Defaults
""""""""

None
