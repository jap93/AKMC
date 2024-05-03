===============
1. Introduction
===============

KMC exploits stochastic algorithms to explore rare events and coarse
grain the time evolution of the model
system. [1]_\ :sup:`,`\  [2]_\ :sup:`,`\  [3]_ In our calculations the
transition states are considered as thermally activated diffusional hops
and are governed by an Arrhenius equation

:math:`r_{D} = D_{0}exp( - \frac{Q}{\text{kT}})` [1]

where *r\ D* is the rate of an event, *D\ 0* is an exponential
pre-factor, *Q* is the activation energy of the hop, *k* is Boltzmann's
constant and *T* the temperature. We assume that the exponential
pre-factor D\ :sub:`0` is 1.0e\ :sup:`13` s\ :sup:`-1` which is the
vibrational frequency. In our simulations an event is either a cation or
an anion diffusional hop and the state of the model system is evolved by
choosing one event stochastically, according to the rate of the events
using the following equation.

:math:`\sum_{i = 1}^{m - 1}{r_{i} \leq \rho_{1}}\sum_{j = 1}^{N}r_{j} < \sum_{k = 1}^{m}r_{k}`
[2]

*where m* is the index of the chosen event and *N* is the total number
of possible events. Summation indices *i*, *j*, and *k* denote the
individual events, thus *r\ i*\ is the rate of the event *i*. *ρ\ i* is
a random number evenly distributed over the range [0, 1). This ensures
that faster events have a greater probability of being chosen than
slower events. The aKMC 500 events within each cycle allowing. Once an
event is chosen, the surface is modified to enact the diffusion event.
The simulation time is then advanced by

:math:`\mathrm{\Delta}t = \  - \frac{\ln\rho_{2}}{\sum_{i = 1}^{N}r_{i}}`
[3]

In equation 3, *t* is the elapsed time and *ρ\ 2*\ is a random number
evenly distributed over the range [0, 1).

The list of events can be pre-determined or calculated on-the-fly. A
pre-determined list requires less computation, but requires a prior
knowledge of the evolution of the simulation cell and/or material
structure and *all* the processes in operation over the entire
simulation. ACDC uses the alternative, yet computationally more
demanding, approach to calculate the activation energies on-the-fly. and
the methodology for calculating the transition states is described in
the next section.

===============================
2. Calculation of saddle points
===============================

It is possible to do this with techniques such as molecular
dynamics, [4]_ the dimer method  [5]_ and activation relaxation
technique (ART). [6]_ At the moment the ART method is the sole technique
coded into ACDC (there are the stubbs for the Dimer method) and have
followed the recipe described in reference [7]_ and briefly describe
the computational method below. ART is member of the minimum mode
following methods that employ an eigenvector-following approach to
efficiently determine transition states. As the transition state is
identified by the lowest eigenvalue-eigenvector pair (λ\ :sub:`1`,
**v**\ :sub:`1`) a Lanczos scheme can be used to iteratively construct a
set of tridiagonal matrices that range from **T**\ :sub:`2` to
**T**\ *j* and whose lowest eigenvalue converges to λ\ :sub:`1` as *j*
increases.[7]_ As many of our simulation cells contain many thousands
of atoms, the Hessian matrix would be prohibitively large to store and
diagonalise. Indeed the tridiagonal form of the matrix can easily be
computed using the equations 1 to 10 in ref [7]_) and diagonalised
(DSTEV from the LAPACK library [ref]). A modified force is then
determined depending whether λ\ :sub:`1` is positive or negative and
used to “push” the system towards a saddle point on the potential energy
surface.

In our calculations we used the FIRE method [8]_ to converge the atomic
positions to the transition state.

===============
3. Installation
===============

The most convenient way to run the program is using empirical potentials. The code can be built in the following way::
cd AccDynamics
cmake -DMLPOT=OFF -DMETPOT=OFF/ON  (when both are OFF the rigid ion model will be used. If METPOT is ON the manybody/metal potentials)
cmake --build . --clean-first

==================
4. Parallelisation
==================

The evaluation of energies and forces are evaluated using OpenMP whilst each transition point search is task farmed using MPI.

=========================
5. How to run the program
=========================


The program requires three input files, *control*, *potentials* and *basis* that contain the
keywords for the functionality of the program, empirical potential parameters and the atomic
positions respectively.

-----------
5.1 control
-----------

The input is broken into sections depending on the functionality. The
main input key words are

+--------------+-------------+------------------------------------------------------+
| **Key word** | **Type**    | **Functionality**                                    |
+==============+=============+======================================================+
| seed         | int         | Seed for the random number. If no number is          |
|              |             | specified the system clock is used. The workgroup    |
|              |             | number is added automatically to seed so that each   |
|              |             | workgroup will follow a different trajectory.        |
+--------------+-------------+------------------------------------------------------+
| parentcores  | int         | The number of mpi threads used by the parent.        |
+--------------+-------------+------------------------------------------------------+
| numgroups    | int         | The number of work groups required. Default is 0.    |
+--------------+-------------+------------------------------------------------------+
| prerelax     |             | Relax the geometry before the KMC starts             |
+--------------+-------------+------------------------------------------------------+
| restart      |             | Restart the calculation with new KMC iteration       |
+--------------+-------------+------------------------------------------------------+
| freeze       | int         | Freeze the following types of atoms. E.g.            |
|              |             |                                                      |
|              |             | freeze 1                                             |
|              |             | Pt                                                   |
|              |             |                                                      |
|              |             | This prevents any movement of Pt type atoms          |
+--------------+-------------+------------------------------------------------------+
| externalfile | string      | The name of any external input file (depricated)     |
+--------------+-------------+------------------------------------------------------+
| rootname     | string      | The working directory (depricated)                   |
+--------------+-------------+------------------------------------------------------+
| close        |             | Finish all input                                     |
+--------------+-------------+------------------------------------------------------+
| geomformat   | string      | The style of geometry input. dlpoly or onetep        |
|              |             | options are available at the moment                  |
+--------------+-------------+------------------------------------------------------+
| kmc{         |             | Starts the processing of the KMC specific key words  |
+--------------+-------------+------------------------------------------------------+
| art{         |             | Starts the processing of the ART specific key words  |
+--------------+-------------+------------------------------------------------------+
| relax{       |             | Starts the processing of the energy minimisation     |
|              |             | specific key words                                   |
+--------------+-------------+------------------------------------------------------+

Description of keywords to control the kinetic Monte Carlo functionality. This input
is started with kmc{ and ended with }. NB the enrgies do not have units but depend on
*energyunit* above. 

+----------------+-------------+------------------------------------------------------+
| **Key word**   | **Type**    | **Functionality**                                    |
+================+=============+======================================================+
| kmcsteps       | int         | The number of kmc steps/cycles                       |
+----------------+-------------+------------------------------------------------------+
| mincap         | double      | The minimum activation energy.                       |
+----------------+-------------+------------------------------------------------------+
| prefactor      | double      | The prefactor used to calculate rates.               |
+----------------+-------------+------------------------------------------------------+
| window         | double      | The window for accepting kmc energies.               |
+----------------+-------------+------------------------------------------------------+
| kmctemperature | double      | The temperature to be used by the kmc simulation     |
+----------------+-------------+------------------------------------------------------+
| kmcevents      | int         | The number of events collected per cycle.            |
+----------------+-------------+------------------------------------------------------+
| basin2delta    | double      | The magnitude of the displacement away from the      |
|                |             | saddle point prior to the attempted relaxation to    |
|                |             | the second basin.                                    |
+----------------+-------------+------------------------------------------------------+
| kmcmethod      | string      | The method used for the saddle search. Only ART is   |
|                |             | available at the moment                              |
+----------------+-------------+------------------------------------------------------+
| kmcbasinradius | double      | The minimum displacement of an atom before it is     |
|                |             | considered to have entered a new basin               |
+----------------+-------------+------------------------------------------------------+
| writeevents    |             | Write all the events out to a file                   |
+----------------+-------------+------------------------------------------------------+
| usegaussian    | double      | The atoms are given a random displacement weighted   |
|                |             | by a Gaussian of the specified width                 |
+----------------+-------------+------------------------------------------------------+

Description of keywords to control the Activation Relaxation Technique. Again units are specified
in the main directives above.

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| numvectors          | int         | The number of Lnanczos vectors used to obtain        |
|                     |             | eigenvalues. Default 20                              |
+---------------------+-------------+------------------------------------------------------+
| maxeigenvalue       | double      | The max eigenvalue. Once an eigenvalue falls below   |
|                     |             | this value the forces parallel to the eigenvalue are |
|                     |             | used.                                                |
+---------------------+-------------+------------------------------------------------------+
| initialdisplacement | double      | The displacement used to activate the ions at the    |
|                     |             | start of the search.                                 |
+---------------------+-------------+------------------------------------------------------+
| eigentolerence      | double      | The tolerence to converge the eigenvalues.           |
+---------------------+-------------+------------------------------------------------------+
| lanczosdisplacement | double      | The displacement of atoms used to calculate the      |
|                     |             | eigenvalues from the tri-diagonal matrix             |
+---------------------+-------------+------------------------------------------------------+
| maxstep             | int         | The number of iterations to calculate the saddle     |
|                     |             | point.                                               |
+---------------------+-------------+------------------------------------------------------+
| minmethod           | string      | The mminimisation technique to find the saddle point |
|                     |             | Only FIRE is available at the moment.                |
+---------------------+-------------+------------------------------------------------------+
| timestep            | double      | The timestep used by FIRE. Typically should be       | 
|                     |             | similar to that used by MD.                          |
+---------------------+-------------+------------------------------------------------------+
| alpha               | double      | The value of alpha used in FIRE                      |
+---------------------+-------------+------------------------------------------------------+
| damp                | double      | damping factor for the parallel forces               |
+---------------------+-------------+------------------------------------------------------+
| forcetol            | double      | The convergence criterion for the minimisation       |
+---------------------+-------------+------------------------------------------------------+

Input keywords for the relaxation

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| debug               |             | increases the amount of information output to files. |
+---------------------+-------------+------------------------------------------------------+
| relaxsteps          | int         | The number of iterations of the minimisaer.          |
+---------------------+-------------+------------------------------------------------------+
| initialdisplacement | double      | The displacement used to activate the ions at the    |
|                     |             | start of the search.                                 |
+---------------------+-------------+------------------------------------------------------+
| maxstep             | double      | The maximum size of the displacement in FIRE.        |
+---------------------+-------------+------------------------------------------------------+
| method              | string      | The mminimisation technique. There is a chiice       |
|                     |             | between FIRE, FIRE2 and external at the moment.      |
|                     |             | (The latter uses the minimisation technique from the |
|                     |             | librray program e.g. onetep.)                        |
+---------------------+-------------+------------------------------------------------------+
| timestep            | double      | The timestep used by FIRE. Typically should be       | 
|                     |             | similar to that used by MD.                          |
+---------------------+-------------+------------------------------------------------------+
| alpha               | double      | The value of alpha used in FIRE                      |
+---------------------+-------------+------------------------------------------------------+
| forcetol            | double      | The convergence criterion for the minimisation       |
+---------------------+-------------+------------------------------------------------------+

--------------
5.1 potentials
--------------

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Type**    | **Functionality**                                    |
+=====================+=============+======================================================+
| cutoff              | double      | The short range cutoff for the potential and Ewald . |
+---------------------+-------------+------------------------------------------------------+
| noimage             |             | The nearest image conventin is used by default. This |
|                     |             | keyword uses a slower multiple image method that     |
|                     |             | is better suited to small simulation cells.          |
+---------------------+-------------+------------------------------------------------------+
| noewald             |             | The Ewald sum is used by default for the two-body    |
|                     |             | potential model. Thus this keyword switches it off.  |
+---------------------+-------------+------------------------------------------------------+
| ewald precision     |             | Controls the accuracy of the Ewald sum. Default      |
|                     |             | value: 1.0e-6                                        |
+---------------------+-------------+------------------------------------------------------+
| species             | int         | species keyword followed by the number of different  |
|                     |             | speccies. Each element type should be input as       |
|                     |             | follows:                                             |
|                     |             | name  mass charge atomic_number                      |
+---------------------+-------------+------------------------------------------------------+

As described in the installation section the program can be compiled with either two-body (including Ewald sum),
many-body (metal potentials) or SchNet potentials. Note all parameters are in electron volts! 
Parameters compatible with the two-body are:

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Parameters**                                                     |
+=====================+=============+======================================================+
| buck                | *A* , :math:`{\alpha}` , *C*                                       |
+---------------------+--------------------------------------------------------------------+
| morse               | *D* , r\ :sub:`eq` , *k*                                           |
+---------------------+--------------------------------------------------------------------+
| ljones              | :math:`{\eta}` , :math:`{\sigma}`                                  |
+---------------------+--------------------------------------------------------------------+
| bhm                 | *A* , :math:`{\alpha}` , *C* , *D*                                 |
+---------------------+--------------------------------------------------------------------+

Here is an example::

   cutoff 8.0
   noimage
   species 2
   Mg 24.0 2.0 12
   O 16.0 -2.0  8
   twobody 2
   buck
   Mg O  1428.5 0.2945  0.00
   buck
   O  O 22764.3 0.1490 27.879
   close

Parameters compatible with the metal potentials are:

+---------------------+-------------+------------------------------------------------------+
| **Key word**        | **Parameters**                                                     |
+=====================+=============+======================================================+
| stch                | :math:`{\eta}` , *a* , *n* , *m* , *c*                             |
+---------------------+--------------------------------------------------------------------+
| gupta               | *A* , r\ :sub:`eq` , *p*, *B*, *q*                                 |
+---------------------+--------------------------------------------------------------------+
| fnsc                | *c0* , *c1* , *c2* , *c* , *A* , *d* , :math:`{\Beta}`             |
+---------------------+--------------------------------------------------------------------+

Here is an example::

   cutoff 6.5
   species 1
   Al  25.0  0.0 13
   manybody 1 ev
   suttonchen
   Al  Al   0.033147    4.05       7.0        6.0         16.399
   close


---------
5.2 basis
---------

The basis format follows that of DL_POLY with a single exception (the number of atoms is input twice). The format is::

   title
   0   3  number_of_atoms  number_of_atoms
   lattice vectors * 3
   name
      x  y  z
   name
      x  y  z

For example::

   CeO2 doped with Gd - bulk
   0    3   2549  2549
   32.7232154959   0.0000000000   0.0000000000
   0.0000000000  32.7232154959   0.0000000000
   0.0000000000   0.0000000000  32.7232154959
   O      1
      1.4168603868  1.3899782202  1.3108695097
   O      2
      4.1014301991  4.0478192159  1.4171530724
   O      3
      3.7163970637  -1.3764087458  4.0580241564

=============
6. References
=============

.. [1]
   D.T. Gillespie, *J. Phys. Chem.*, 1997, **81**, 2340-2361.

.. [2]
   A.F. Voter, *Phys. Rev. B*, 1986, **34**, 6819-6829.

.. [3]
   C.C. Battaile, *Comput. Methods Appl. Mech. Engrg.*, 2008, **197**,
   3386-3398.

.. [4]
   D. Frenkel and B. Smit, *Understanding Molecular Simulation: From
   Algorithms to Applications*, 2002, Academic Press.

.. [5]
   G. Henkelman and H. Jónsson, *J. Chem. Phys.*, 1999, **111**,
   7010-7020.

.. [6]
   G.T. Barkema and N. Mousseau, *Computational Materials Science*,
   2001, **20**, 285-292.

.. [7]
   R.A. Olsen, G.J. Kroes, and G. Henkelman, *The Journal of Chemical Physics*, 2004,  **121(20)**, 9776–9792.

.. [8]
   E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, P. Gumbsch, *Phys. Rev. Lett.*, 2006, **97** , Article 170201.
