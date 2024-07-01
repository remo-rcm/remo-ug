.. _`chap:remo`:

==============
The REMO Model
==============

Throughout this study new methodical developments will be implemented
and validated with the regional hydrostatic climate model REMO REMO. In
this chapter a brief overview of the REMO model in general and its
dynamical core in particular is given.

.. _`sec:remo_overview`:

Overview
--------

REMO :raw-latex:`\parencite{Jacob2001,Jacob2001b}` is a
three-dimensional hydrostatic primitive equation regional climate model
originally based on the German Weather Service “Europa-Modell”
:raw-latex:`\parencite{Majewski1991}`. To allow application to climate
projections the physics package of the Max Planck Institute for
Meteorology global climate model ECHAM4
:raw-latex:`\parencite{Roeckner1996}` was implemented into REMO and has
since been continuously developed and expanded. More recently the model
has also been updated with a non-hydrostatic extension
:raw-latex:`\parencite{Goettel2008}`.

.. _`subsec:remo_physics`:

Physics
~~~~~~~

The present work will exclusively deal with the dynamical component of
the model, i.e., the solution of the primitive equations in Euler-form
and more specifically the computation of the pressure gradient force.
Giving even a cursory overview of the numerous physical
parameterizations in REMO, i.e., those parts of the model that describe
diabatic processes such as, e.g., freezing and melting, microphysical
interactions, radiative transfer and other components of the climate not
covered by the primitive equations, is hence beyond the scope of this
work. Instead some of the main parameterization schemes are listed and
appropriate references are pointed out for detailed descriptions. The
main characteristics of the physical parameterizations are given by a
radiation scheme based on :raw-latex:`\textcite{Morcrette1986}` and
:raw-latex:`\textcite{Giorgetta1995}`. The treatment of stratiform
clouds is governed by schemes based on
:raw-latex:`\textcite{Sundqvist1978}` with modifications based on
:raw-latex:`\textcite{Roeckner1996}`, where also a lot of in-depth
information on several other parameterizations in REMO can be found.
Cumulus convection is represented through a mass flux scheme after
:raw-latex:`\textcite{Tiedtke1989}` and
:raw-latex:`\textcite{Nordeng1994}`. The turbulent surface fluxes are
computed according to Monin-Obukhov similarity theory
:raw-latex:`\parencite{Louis1979}` with vertical diffusion based on the
turbulent kinetic energy. Soil processes are modeled with a 5-layer
diffusion scheme for heat transfer and a bucket approach for soil
moisture, interception by vegetation and snow with modifications that
allow the freezing and thawing of soil water
:raw-latex:`\parencite{Semmler2002}`. The runoff scheme is based on
:raw-latex:`\textcite{Duemenil1992}`.

.. _`subsec:remo_dynamics`:

Dynamics
~~~~~~~~

The most important aspects of REMO’s dynamical core will be covered in
the following sections in some detail, but here the main features as
relevant for the developments in later chapters are summarized. The
basic model formulation is based on the hydrostatic primitive equation
system :raw-latex:`\parencite[e.g.,][]{White2005}`), i.e., a rendition
of the Euler-equations of motion under the assumption of a
hydrostatically balanced atmosphere. In such cases the third equation of
motion reduces to a diagnostic relationship for the geopotential height,
significantly reducing computational costs in the process, but also
limiting the model to coarser resolutions beyond about ten kilometers.

The basic equations are discretized using a centered finite-difference
approach on a regular longitude-latitude grid in the horizontal and a
hybrid terrain-following grid in the vertical. For the time
discretization the leapfrog scheme is employed, necessitating the use of
a filter to maintain stability of the computation. In REMO the
Asselin-filter is used for this purpose. While centered differences in
general would allow for second order approximations, in practice the use
of the Asselin-filter limits the order of approximation somewhat below
this mark :raw-latex:`\parencite[see e.g.,][]{Duran2010}`.

Stability requirements also motivate the addition of diffusion terms to
the equations in order to prevent the accumulation of numerical noise at
the shortest wavelengths resolved by the model. To lessen the
constraints given by the CFL number
:raw-latex:`\parencite[e.g.,][]{Duran2010}`, some of the model also
features implicit formulations. On the one hand the vertical advection
is handled by the classical Crank-Nicolson scheme. On the other hand
processes prone to generating gravity waves are dealt with applying the
scheme of :raw-latex:`\textcite{Simmons1981}` but can also be approached
explicitly instead, depending on model configuration. Implicit
formulations enable the use of larger timesteps, generally outweighing
the drawbacks of having to solve additional equations in each of them.

Furthermore, since REMO is a regional model boundary conditions need to
be provided in addition to an initial state, generally created from
global model output or reanalysis products. To smoothly impose the
boundary values on the solution a lateral relaxation scheme after
:raw-latex:`\textcite{Davies1976}` is applied.

.. _`sec:remo_prelim`:

Coordinate Systems
------------------

In this section the coordinate systems used in REMO are described as a
prerequisite to formulating the equations of motions later on.

.. _`subsec:horiz_coord_system`:

Horizontal Coordinate System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REMO employs a spherical coordinate system in the horizontal directions.
It is derived by intricate rotation of the standard geographical
coordinate system. The rotation is constructed such that the rotated
equator passes through the center of the domain. Near the poles
meridional convergence can result in very small computational cells,
effectively dominating the choice of timestep for the whole domain due
to the CFL-condition :raw-latex:`\parencite[e.g.,][]{Duran2010}`.
Because REMO is a regional model the rotation in general guarantees that
the rotated domain is far away from the poles, as illustrated by
`1.1 <#fig:rotation>`__. Consequently, the influence of meridional
convergence on the choice of viable timesteps is strongly limited. A
detailed account of the coordinate transformation, its inverse and the
conversion of wind components between the rotated and the geographical
system can be found in the literature :raw-latex:`\parencite{emdm}`.
Note however that the resulting equation systems – before and after
transformation – are exactly the same with only one exception: the
Coriolis factor :math:`f = 2\Omega\sin{\varphi}`. In the rotated system
one instead has:

.. math::

   \label{eq:coriolis}
       f = 2\Omega\left(\sin\varphi\sin\varphi^N + \cos\varphi\cos\varphi^N\cos\left(\lambda-\lambda^N\right)\right)

where :math:`\left(\lambda, \varphi\right)` are the rotated longitude
and latitude respectively and :math:`\left(\lambda^N, \varphi^N\right)`
are the coordinates of the geographical north pole in the rotated
system.

.. container:: float
   :name: fig:rotation

   .. container:: center

      |image1|

.. _`subsec:vert_coord_system`:

Vertical Coordinate System
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the vertical direction the REMO model uses a terrain-following
coordinate system. The classical example of such a coordinate for
meteorological use is the :math:`\sigma`-coordinate developed by
:raw-latex:`\textcite{Phillips1957}`. The idea is to express the model
height not in terms of, e.g., geometric or geopotential height, or in
terms of the atmospheric pressure but as the ratio of hydrostatic
pressure to hydrostatic surface pressure:

.. math:: \sigma = \frac{p}{p_s}

Clearly at surface level one has :math:`p = p_s` and hence
:math:`\sigma = 1`. The lowest coordinate surface in the
:math:`\sigma`-system coincides with the bottom topography. On the other
hand at the top of the atmosphere one has :math:`p = 0` and hence also
:math:`\sigma = 0`. Since the coordinate is built from the hydrostatic
pressure – which of course is a monotone function of height – the result
is also a monotone function, allowing its use as a vertical coordinate
in the equations of motion. The main advantage is that the lower
vertical boundary condition simplifies substantially, avoiding the use
of, e.g., uncentred differences in the presence of sloping topography.

There are however also a number of tradeoffs: meteorological
observations often are given on surfaces of constant pressure. Such data
has to be interpolated to :math:`\sigma`-surfaces to be used as driving
data in terrain-following models (e.g.,
:raw-latex:`\cite{Sundqvist1976}`). Additionally while the sloping of
coordinate surfaces is useful at the lower boundary, it serves no
purpose in the upper atmosphere. Here typical flow regimes tend to
follow isobaric surfaces that would be most naturally expressed in a
pressure coordinate. Consequently, the representation along
:math:`\sigma`-surfaces is often numerically more difficult and can
produce noise around steeply sloped coordinate surfaces even in the
upper atmosphere.

Moreover, the representation of the horizontal pressure gradient force
is also more complicated in this and related coordinates
:raw-latex:`\parencite[e.g.,][]{Sundqvist1978}`. To reduce these
disadvantages REMO instead uses a hybrid coordinate after
:raw-latex:`\textcite{Simmons1981}`, which will be referred to as the
:math:`\eta`-coordinate during this work. The idea is to use a pressure
coordinate in a large part of the upper atmosphere thereby facilitating
a natural representation of the flow. At the lower boundary a
:math:`\sigma`-system is used instead to retain the advantages of
terrain following coordinates in regards to the lower boundary
condition. Between the two regimes the coordinate is linearly
interpolated to allow a suitable transition.

This can be written

.. math::

   \label{eq:eta}
       \eta =
       \begin{cases}
           \frac{p}{p_r}                                                     & \text{for}\ 0\leq p\leq p_t   \\
           \frac{p-p_t}{p_s-p_t} + \frac{p_t}{p_r}\cdot\frac{p_s-p}{p_s-p_t} & \text{for}\ p_t\leq p\leq p_s
       \end{cases}

with :math:`p_r` a constant reference pressure and :math:`p_t` a
pressure threshold above which :math:`\eta` is identical to a pressure
coordinate. In REMO :math:`p_r = \SI{1013.25}{hPa}` is used for the
reference pressure. The pressure threshold varies between simulations
but is generally chosen roughly around :math:`\SI{220}{hPa}`. From
`[eq:eta] <#eq:eta>`__ it can be readily seen that one has
:math:`\eta = 1` for :math:`p = p_s`, i.e., the lowest coordinate
surface exactly follows the given orography. As will be seen later this
is enough to ensure the benefits of a simplified lower boundary
condition. At :math:`p = p_t` the two terms in `[eq:eta] <#eq:eta>`__
are exactly identical and the transition between the terrain-following
and the pressure regime is continuous. Clearly, the influence of
orographic variation on the coordinate surfaces is decaying linearly
towards the threshold :math:`p_t`. For any pressure :math:`p\leq p_t`
the coordinate surfaces are indeed completely independent of the given
surface elevation. In this way detrimental numerical effects of the
sloping :math:`\eta`-surfaces as discussed above are effectively
limited. Furthermore, the pressure surfaces in the upper atmosphere are
well suited to describe flow in the free atmosphere.
`1.4 <#fig:lorenz>`__ illustrates these typical features of the
coordinate surfaces in the :math:`\eta`-system.

According to `[eq:eta] <#eq:eta>`__ :math:`\eta` is a function of
:math:`p` and :math:`p_s`, but likewise one can express :math:`p` as a
function of :math:`\eta` and :math:`p_s`. With the definitions

.. math::

   \label{eq:ak}
       A(\eta) =
       \begin{cases}
           p_r\eta                                    & \text{for}\ 0\leq\eta\leq\eta_t  \\
           \frac{p_r p_t}{p_r-p_t}\left(1-\eta\right) & \text{for}\ \eta_t\leq\eta\leq 1
       \end{cases}

.. math::

   \label{eq:bk}
       B(\eta) =
       \begin{cases}
           0                           & \text{for}\ 0\leq\eta\leq\eta_t  \\
           \frac{p_r\eta-p_t}{p_r-p_t} & \text{for}\ \eta_t\leq\eta\leq 1
       \end{cases}

where :math:`\eta_t` is the :math:`\eta`-value corresponding to the
pressure threshold :math:`p_t`, i.e., :math:`\eta_t = \frac{p_t}{p_r}`,
we can write:

.. math::

   \label{eq:akbkps}
       p = A(\eta) + B(\eta)\cdot p_s

This means that in the :math:`\eta`-system the atmospheric pressure
takes the form of a linear function. As will be seen later equation
`[eq:akbkps] <#eq:akbkps>`__ can be used as a diagnostic equation to
recover the pressure from the given atmospheric parameters
:math:`A(\eta), B(\eta)` and the hydrostatic surface pressure
:math:`p_s`.

Continuous Model Equations
--------------------------

In this section the basic continuous equations of motion, used to
advance the meteorological variables of interest in time, are given.
This exposition closely follows :raw-latex:`\textcite{emdm}` where more
in-depth information about the dynamical core can be found.

Prognostic Equations
~~~~~~~~~~~~~~~~~~~~

Here the main equations for each of the prognostic variables are
formulated: hydrostatic surface pressure, horizontal wind components,
temperature as well as specific humidity, cloud water and cloud ice.

Equation of Continuity
^^^^^^^^^^^^^^^^^^^^^^

Technically the equation of continuity is not used as a prognostic
equation for the REMO model, but the equations for the surface pressure
and the vertical velocity in the :math:`\eta`-system both are derived
from it. For completeness it is hence given here as

.. math::

   \label{eq:cont}
       \dfdt{}\left(\dfdeta{p}\right) + \acosphi\left[\dfdlam{}\left(u\dfdeta{p}\right) + \dfdphi{}\left(v\cos{\varphi}\dfdeta{p}\right)\right] + \dfdeta{}\left(\dot{\eta}\dfdeta{p}\right) = 0

where :math:`p` is the hydrostatic pressure, :math:`a = \SI{6371229}{m}`
the radius of earth, :math:`\varphi` and :math:`\lambda` the rotated
latitude and longitude respectively, :math:`u, v` the zonal and
meridional wind components, and :math:`\dot{\eta} = \frac{D\eta}{Dt}`
the vertical velocity in the :math:`\eta`-system. As pointed out earlier
one of the main advantages of a terrain-following coordinate system is
the simple formulation of vertical and in particular the lower boundary
condition. Assuming that no mass transport happens through the upper or
lower boundary, i.e., the interfaces to outer space and inner earth, the
boundary conditions can simply be given as:

.. math::

   \begin{aligned}
       \label{eq:vbc}
       \begin{split}
       \etadot &= 0\ \text{at the upper boundary with}\ \eta = 0 \\ 
       \etadot &= 0\ \text{at the lower boundary with}\ \eta = 1
       \end{split}
   \end{aligned}

.. _`ssubsec:remo_ps`:

Surface Pressure
^^^^^^^^^^^^^^^^

Integrating `[eq:cont] <#eq:cont>`__ from top to bottom of the
atmosphere and employing the vertical boundary conditions
`[eq:vbc] <#eq:vbc>`__ yields a prognostic equation for the hydrostatic
surface pressure :math:`p_s`. The rate of change of :math:`p_s` is then
given by:

.. math::

   \label{eq:ps}
       \frac{\partial p_s}{\partial t}
       + \frac{1}{a\cos\varphi}\intdeta{0}{1}{\frac{\partial}{\partial\lambda}\left(u\frac{\partial p}{\partial\eta}\right)
           + \frac{\partial}{\partial\varphi}\left(v\cos\varphi\frac{\partial p}{\partial\eta}\right)}
       = 0

.. _`ssubsec:remo_uv`:

Wind Components
^^^^^^^^^^^^^^^

The evolution of the wind component in zonal direction :math:`u` is
described by the following equation:

.. math::

   \label{eq:u}
       \frac{\partial u}{\partial t} -
       \frac{1}{\cos\varphi}Q\frac{\partial p}{\partial\eta}v\cos\varphi
       + \frac{1}{a\cos\varphi}\frac{\partial}{\partial\lambda}\left(\Phi+K\right)
       + \frac{RT_v}{a\cos\varphi}\frac{\partial\ln p}{\partial\lambda}
       + \dot\eta\frac{\partial u}{\partial\eta}
       - D_u
       = F_u

Here :math:`Q` denotes the absolute potential vorticity, :math:`p` the
hydrostatic pressure, :math:`\Phi` the geopotential, :math:`T_v` the
virtual temperature, :math:`R = \SI{287.05}{J/(kg.K)}` the gas constant
for an ideal gas, :math:`\etadot` the vertical velocity and :math:`D_u`
the horizontal diffusion of :math:`u`. Likewise for the the wind
component in meridional direction :math:`v` one has

.. math::

   \label{eq:v}
       \frac{\partial v}{\partial t} -
       Q\frac{\partial p}{\partial\eta}u
       + \frac{1}{a}\frac{\partial}{\partial\varphi}\left(\Phi+K\right)
       + \frac{RT_v}{a}\frac{\partial\ln p}{\partial\varphi}
       + \dot\eta\frac{\partial v}{\partial\eta}
       - D_v
       = F_v

where :math:`D_v` is the horizontal diffusion of :math:`v`. The terms
:math:`F_u` and :math:`F_v` represent diabatic and subscale, i.e.,
parameterized processes such as convection or turbulent boundary layer
interactions. Since physical parameterizations are not the topic of this
work and for reasons of brevity the reader is referred to the various
sources cited in `1.1.1 <#subsec:remo_physics>`__ for detailed
treatments. This also applies to related terms appearing in the
remaining equations.

.. _`ssubsec:remo_T`:

Temperature
^^^^^^^^^^^

The thermodynamic equation is given by

.. math::

   \label{eq:T}
       \dfdt{T} + \acosphi\left(u\dfdlam{T} + v\cosphi\dfdphi{T}\right) + \etadot\dfdeta{T} - D_T = \frac{\alpha\omega}{c_p} + F_T

where :math:`T` is the absolute temperature, :math:`\alpha` the specific
volume, :math:`\omega = \frac{dp}{dt}` the vertical velocity in terms of
the hydrostatic pressure :math:`p`, :math:`D_T` the horizontal diffusion
of temperature and :math:`F_T` changes in temperature due to subscale
processes.

.. _`ssubsec:remo_Q`:

Specific Humidity, Cloud Water and Cloud Ice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The moist components of air are realized as passive tracers

.. math::

   \label{eq:qd}
       \dfdt{q_d} + \acosphi\left(u\dfdlam{q_d} + v\cosphi\dfdphi{q_d}\right) + \etadot\dfdeta{q_d} - D_{q_d} = F_{q_d}

.. math::

   \label{eq:qw}
       \dfdt{q_w} + \acosphi\left(u\dfdlam{q_w} + v\cosphi\dfdphi{q_w}\right) + \etadot\dfdeta{q_w} - D_{q_w} = F_{q_w}

.. math::

   \label{eq:qi}
       \dfdt{q_i} + \acosphi\left(u\dfdlam{q_i} + v\cosphi\dfdphi{q_i}\right) + \etadot\dfdeta{q_i} - D_{q_i} = F_{q_i}

with specific humidity :math:`q_d`, cloud water :math:`q_w` and cloud
ice :math:`q_i` and their horizontal diffusion terms
:math:`D_{q_d}, D_{q_w}, D_{q_i}` respectively. Feedback mechanisms are
included by using the virtual temperature :math:`T_v` in
`[eq:u] <#eq:u>`__ and `[eq:v] <#eq:v>`__ instead of the absolute
temperature :math:`T`. As before :math:`F_{q_d}, F_{q_w}` and
:math:`F_{q_i}` represent subscale interactions including for example
phase conversions.

Diagnostic Equations
~~~~~~~~~~~~~~~~~~~~

Equations `[eq:ps] <#eq:ps>`__ to `[eq:qi] <#eq:qi>`__ constitute a
closed system that can in principle be solved for the prognostic
variables :math:`p_s, u, v, T, q_d, q_w` and :math:`q_i`. First however,
the auxiliary quantities – such as for instance the geopotential
:math:`\Phi` or the absolute potential vorticity :math:`Q` – appearing
in these equations have to be expressed in terms of the prognostic
variables. Therefore, in this section the definitions of these
diagnostic quantities will be introduced.

Pressure
^^^^^^^^

Equation `[eq:ps] <#eq:ps>`__ references the pressure :math:`p` at
arbitrary levels. In hybrid pressure based coordinates the pressure at a
given level is a linear function of the surface pressure :math:`p_s` and
the vertical coordinate :math:`\eta`. With the definitions of
`1.2.2 <#subsec:vert_coord_system>`__ one can write:

.. math:: p = A(\eta) + B(\eta)\cdot p_s

Potential Absolute Vorticity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The potential absolute vorticity :math:`Q` is needed to evaluate the
momentum equations `[eq:u] <#eq:u>`__ and `[eq:v] <#eq:v>`__. It is
given by

.. math::

   \label{eq:pav}
       Q = \left(\dfdeta{p}\right)^{-1}\left(f + \acosphi\left(\dfdlam{v}-\dfdphi{u\cos\varphi}\right)\right)

where :math:`f` is the Coriolis force. According to
`[eq:coriolis] <#eq:coriolis>`__ this term can be written

.. math:: f = 2\Omega\left(\sin\varphi\sin\varphi^N + \cos\varphi\cos\varphi^N\cos\left(\lambda-\lambda^N\right)\right)

with

.. math:: \Omega = \SI{7.29211e-5}{s^{-1}}

the angular velocity.

Geopotential
^^^^^^^^^^^^

In a hydrostatic model such as REMO the geopotential is characterized by
the balance of a gravity and a buoyancy term. Essentially this balance
is the residual of the third equation of motion under hydrostatic
conditions and given by

.. math::

   \label{eq:hydro_eta}
       -RT_v\dfdeta{\ln{p}} = \dfdeta{\Phi}

and consequently the geopotential can be recovered by vertical
integration from the bottom of the atmosphere to a given height
:math:`\eta`:

.. math::

   \label{eq:geop}
       \Phi = \Phi_s - R\intdeta{1}{\eta}{T_v\dfdeta{\ln{p}}}

The term :math:`\Phi_s = g\cdot z_s` denotes the surface geopotential
with the height above mean sea level :math:`z_s`.

Kinetic Energy
^^^^^^^^^^^^^^

Term :math:`K` in the momentum equations `[eq:u] <#eq:u>`__ and
`[eq:v] <#eq:v>`__ represents the kinetic energy per unit mass. It is
given by

.. math::

   \label{eq:kin}
       K = \frac{1}{2}\left(u^2 + v^2\right)

Virtual Temperature
^^^^^^^^^^^^^^^^^^^

The virtual temperature :math:`T_v` can be expressed as a function of
the absolute temperature :math:`T` and the moist components of air
:math:`q_d, q_w` and :math:`q_i` as

.. math::

   \label{eq:vtemp}
       T_v = T\left(1 + \left(\frac{R_D}{R} - 1\right)q_d - \left(q_w + q_i\right)\right)

where :math:`R_d = \SI{461.51}{J/(kg.K)}` is the gas constant for water
vapor.

Vertical Velocity in the :math:`\eta`-System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A diagnostic equation for the vertical velocity :math:`\etadot`, i.e.,
the change of the vertical coordinate :math:`\eta` following an air
parcel along its trajectory, can be derived from
`[eq:cont] <#eq:cont>`__. Integrating from top of the atmosphere to a
given height :math:`\eta` yields

.. math::

   \label{eq:etadot}
       \etastar = -\left(\frac{\partial p}{\partial p_s}\right)\dfdt{p_s} - \acosphi\intdeta{0}{\eta}{\frac{\partial}{\partial\lambda}\left(u\frac{\partial p}{\partial\eta}\right)
           + \frac{\partial}{\partial\varphi}\left(v\cos\varphi\frac{\partial p}{\partial\eta}\right)}

with the auxiliary quantity :math:`\etastar = \etadot\dfdeta{p}`.

Specific Volume
^^^^^^^^^^^^^^^

The specific volume :math:`\alpha` can be derived from the equation of
state which can be written as:

.. math::

   \label{eq:alpha}
       \alpha = \frac{R T_v}{p}

Note that for a dry atmosphere `[eq:alpha] <#eq:alpha>`__ simply reduces
to the ideal gas equation.

Vertical Velocity in the :math:`p`-System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :math:`\frac{\alpha\omega}{c_p}` term appearing in the thermodynamic
equation `[eq:T] <#eq:T>`__ is related to the conversion of potential
and kinetic energy. A diagnostic expression for :math:`\omega` can be
found by direct differentiation of :math:`\omega = \frac{Dp}{Dt}` and
substituting equation `[eq:etadot] <#eq:etadot>`__. This yields:

.. math::

   \begin{split}
       \label{eq:omega}
       \frac{\omega}{p} = -\frac{1}{p}\left[\acosphi\intdeta{0}{\eta}{\frac{\partial}{\partial\lambda}\left(u\frac{\partial p}{\partial\eta}\right)
               + \frac{\partial}{\partial\varphi}\left(v\cos\varphi\frac{\partial p}{\partial\eta}\right)}\right] + \\
               \acosphi\left(u\dfdlam{\ln{p}} + v\cosphi\dfdphi{\ln{p}}\right)
   \end{split}

Initial- and Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the definitions of this section equations
`[eq:ps] <#eq:ps>`__–`[eq:qi] <#eq:qi>`__ constitute a closed system
that can in principle be solved for the prognostic variables
:math:`p_s, u, v, T, q_i, q_d, q_w`. This depends of course on the
provision of suitable initial and boundary conditions. For the initial
conditions the values of the prognostic variables are prescribed at the
initial time and for the whole computational domain. In practice
simulations of historical climate or reanalysis data are used to acquire
appropriate initial values. For the vertical boundary an appropriate
condition has already been established with `[eq:vbc] <#eq:vbc>`__. In
case of a global model no further boundary conditions would be required.
REMO however is a regional model and therefore also requires suitable
conditions at the lateral boundaries of the domain. Again such values
can be acquired from driving model simulations, both global and regional
ones, or reanalysis data.

However – unlike for the vertical boundary – there is no obvious way to
impose the resulting values for the prognostic variables near the
lateral boundary. In REMO this is resolved by gradually relaxing the
prognostic variables towards a specified reference state at the boundary
according to :raw-latex:`\textcite{Davies1976}`. For any prognostic
variable :math:`\psi` with given reference state :math:`\psi_R` at the
lateral boundary, an additional relaxation term of the form
:math:`\mu_R\left(\psi-\psi_R\right)` results on the right hand side of
the respective prognostic equation
`[eq:ps] <#eq:ps>`__–`[eq:qi] <#eq:qi>`__. The relaxation factor
:math:`\mu_R` is chosen in such a way that one has :math:`\mu_R = 1` at
the lateral boundary and a rapid but smooth decay towards zero outside
of the boundary.

Consequently, in practice :math:`\mu_R` is essentially zero everywhere
but in a boundary zone of generally less than a dozen grid cells. In
this way the driving values are imposed only at the boundary. At the
same time the boundary zone allows for a smooth transition between the
inner domain characterized by regional dynamics and the lateral boundary
dominated by the forcing patterns. This strategy is effective but also
likely to cause numerical noise near the lateral boundaries. In practice
simulation results within the boundary zone are therefore treated
cautiously and often excluded from analysis for increased robustness.

Discretization
--------------

In the last section the continuous equations of motions in the hybrid
terrain-following :math:`\eta`-coordinate were given and appropriate
initial and boundary conditions were specified. Mathematically this
should ensure the existence and uniqueness of a solution for which the
given equation system `[eq:ps] <#eq:ps>`__–`[eq:qi] <#eq:qi>`__ can then
be solved. However, in general complex partial differential equations
can not be solved explicitly, i.e., in closed form
:raw-latex:`\parencite[e.g.,][]{Duran2010}`. Instead solutions to such
systems are typically approximated numerically, e.g., on high
performance computing systems owing to the complexity of the task. Any
computer can by its very nature only represent a finite number of
distinct states, e.g., due to finite memory and a number of other
limitations. On the other hand the system
`[eq:ps] <#eq:ps>`__–`[eq:qi] <#eq:qi>`__ – formulated in a mathematical
continuum – can assume an infinite number of states.

In other words: computer systems are in general unable to solve such
systems directly. Instead the equations first have to be rendered into a
discrete form, i.e., one that requires only a finite number of unknowns
to represent the state of the system across the domain of interest.
Naturally, this representation will – except for very simple states – be
only an approximation to the true solution of the underlying equation
system. The more discrete points are used to represent the system, the
better this approximation will be, as illustrated in
`1.2 <#fig:discretization>`__. Each gridbox shown is represented in the
discrete model by a single value for each of the prognostic variables.
With more gridboxes per area local structures including topographic
features can be resolved much better. In the following equations
`[eq:ps] <#eq:ps>`__–`[eq:qi] <#eq:qi>`__ will be cast into a discrete
form that facilitates the approximation of solutions on computer
systems.

.. container:: float
   :name: fig:discretization

   .. container:: center

      |image2|

Grid Structure and Discrete Operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, an overview of the discrete structure employed in the REMO model
is given. As mentioned before the idea is to approximate a mathematical
continuum by a finite number of data points distributed within the
region of interest.

Horizontal Grid
^^^^^^^^^^^^^^^

Horizontally REMO employs the so called Arakawa-C grid
:raw-latex:`\parencite[see][]{Arakawa1977}`. The main feature of this
grid is that not all prognostic variables are distributed along the same
discrete points, i.e., the grid is staggered. The main discretization
points are labeled with full indices :math:`(i, j)` as seen in
`1.3 <#fig:arakawac>`__.

.. container:: float
   :name: fig:arakawac

   .. container:: center

      |image3|

At these points (shown as diamond nodes), often referred to as mass
points, all the prognostic variables except for the horizontal wind
components are located. The mass points are placed equidistantly within
the rotated coordinate system with increments of :math:`\Delta\lambda`
and :math:`\Delta\varphi`, where in REMO generally
:math:`\Delta\lambda = \Delta\varphi`. If one thinks of mass points as
the centers of computational cells (dashed lines in
`1.3 <#fig:arakawac>`__) the horizontal wind components :math:`u` and
:math:`v` are placed at the interfaces of those cells. That is, the
:math:`u`-component is displaced to the right by
:math:`\frac{\Delta\lambda}{2}` and the :math:`v`-component upwards by
:math:`\frac{\Delta\varphi}{2}`. Consequently, u-points labeled with
indices :math:`(\iph,j)` and :math:`(i,\jph)` for v-points respectively.
The vorticity :math:`\zeta` is displaced in both horizontal directions
and labeled :math:`(\iph, \jph)`. Note that despite being equidistantly
chosen within the computational space (i.e., the rotated
longitude-latitude system) in physical space (i.e., the geographical
system) the effective horizontal resolution differs between the cells.
This is most notably seen towards the poles due to meridional
convergence.

Vertical Grid
^^^^^^^^^^^^^

Vertically the model is discretized with the so called Lorenz grid.
Again the main feature of this grid is the staggering of its variables
as illustrated in `1.4 <#fig:lorenz>`__. The solid lines are referred to
as full layers and here all prognostic variables are located except for
the (surface) pressure. That (as well as the diagnostic quantities
:math:`\Phi` and :math:`\etadot`) is instead located on half layers
(dashed lines in `1.4 <#fig:lorenz>`__), halfway between the full
layers.

This construction is somewhat similar to the horizontal placement of
variables, but note one of the main differences. In the horizontal
arrangement the discretization points are distributed equidistantly at
least within the rotated coordinate system. In the vertical grid this is
generally not the case. Instead the vertical resolution in the
:math:`\eta`-system can differ substantially between different
computational cells. Essentially the :math:`\eta`-points at which the
variables are to be evaluated are chosen freely and do not change over
the course of a simulation. Note however that again a distinction
between the computational space and the physical space has to be made.
Unlike the horizontal resolution which in physical space varies between
cells but is constant in time, the vertical resolution in physical space
generally changes in every timestep depending on the value of the
surface pressure.

The initial resolution (i.e., based on a reference pressure) in practice
is chosen such that the planetary boundary layer is well-resolved. With
increasing height the resolution becomes much coarser. For instance the
lowest layer in standard setups often is only about 30 m thick, where
the highest layers can reach thicknesses of 3 km and more. For the
labeling of grid indices in the vertical :math:`k` is used for full and
:math:`\kph` for half layers. For instance a discrete temperature value
on the grid is written :math:`T_{ijk}`. Note however that in most cases
one or several of the indices can be inferred from context. Therefore
the same value will often be denoted as :math:`T_{ij}` or :math:`T_k`
instead to improve the readability of formulas.

.. container:: float
   :name: fig:lorenz

   .. container:: center

      |image4|

Finite Difference and Averaging Operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To construct a closed discrete equation system that allows computing the
prognostic variables at the discretization points, the continuous
equations have to be rendered into a discrete form. The main difficulty
here is the approximation of the horizontal and vertical derivatives. In
REMO finite-difference approximations, more specifically so called
centered-differences, are used to achieve this. The main advantage of
centered-differences is that they in principle allow for second order
approximations of the derivatives, i.e., with increasing resolution
:math:`\Delta\lambda, \Delta\varphi` the discretization error will tend
to zero quadratically. At the same time, centered-differences are
comparatively easy to implement in contrast to approximations that allow
even higher orders of approximation. As discussed above the staggered
nature of the computational grid facilitates the use of centered
differences, making them the natural choice for REMO.

However, the use of staggered grids also requires in some places the
transfer of quantities at mass points to :math:`u`- or :math:`v`-points
and vice versa. The same applies to the vertical grid points. To achieve
this a simple averaging approach is used which also retains the order of
approximation.

The difference operators are defined as follows:

.. math::

   \begin{aligned}
       \label{eq:deltaphilam}
       \begin{split}
           \delta_\lambda\psi &= \frac{\psi_{i+1j}-\psi_{ij}}{\Delta\lambda} \\
           \delta_\varphi\psi &= \frac{\psi_{ij+1}-\psi_{ij}}{\Delta\varphi}
       \end{split}
   \end{aligned}

where :math:`\psi` is any quantity located at mass-points, such as
virtual temperature :math:`T_v` or surface pressure :math:`p_s`. Note
that according to equation `[eq:deltaphilam] <#eq:deltaphilam>`__ the
discrete differences of such quantities in :math:`\lambda`-direction are
defined at :math:`u`-points, and the ones in :math:`\varphi`-direction
at :math:`v`-points respectively. At these points the operators in
`[eq:deltaphilam] <#eq:deltaphilam>`__ constitute second-order accurate
centered-differences. For quantities located on :math:`u`- or
:math:`v`-points one likewise defines

.. math::

   \begin{aligned}
       \label{eq:deltaphilamuv}
       \begin{split}
           \delta_\lambda\psi &= \frac{\psi_{\iph j}-\psi_{\imh j}}{\Delta\lambda} \\
           \delta_\varphi\psi &= \frac{\psi_{i\jph}-\psi_{i\jmh}}{\Delta\varphi}
       \end{split}
   \end{aligned}

and note that the operators `[eq:deltaphilamuv] <#eq:deltaphilamuv>`__
are defined at mass-points, i.e., at grid cells labeled with full
indices :math:`(i, j)`.

It is in some cases required to evaluate quantities naturally situated
at :math:`u-` or :math:`v-`\ points at mass points and vice versa. The
following average operators allow for a simple but effective conversion:

.. math::

   \begin{aligned}
       \label{eq:avgmass}
       \begin{split}
           \overline{\psi}^\lambda = \frac{\psi_{i+1 j} + \psi_{i j}}{2} \\
           \overline{\psi}^\varphi = \frac{\psi_{i j+1} + \psi_{i j}}{2}
       \end{split}
   \end{aligned}

Analogously, for those quantities located on u- or v-points one instead
has:

.. math::

   \begin{aligned}
       \label{eq:avguv}
       \begin{split}
           \overline{\psi}^\lambda = \frac{\psi_{\iph j} + \psi_{\imh j}}{2} \\
           \overline{\psi}^\varphi = \frac{\psi_{i\jph} + \psi_{i\jmh}}{2}
       \end{split}
   \end{aligned}

Vertically one has to convert quantities located on half-layers to
full-layers occasionally and this is again achieved by means of vertical
averaging:

.. math::

   \begin{aligned}
       \label{eq:avgvert}
       \begin{split}
           \overline{\psi}^\eta &= \frac{\psi_{\kph}+\psi_{\kmh}}{2}
       \end{split}
   \end{aligned}

Note that generally indices are omitted wherever they can be inferred
from context. For instance, in `[eq:avgvert] <#eq:avgvert>`__ all
references to the horizontal indices :math:`i` and :math:`j` are omitted
as they are not relevant for the vertical averaging.

Leapfrog Time Integration
^^^^^^^^^^^^^^^^^^^^^^^^^

With the previously defined horizontal and vertical difference operators
the spatial derivatives can be approximated. However, time derivatives
also need to be approximated by suitable finite-differences. In REMO the
leapfrog scheme is used for this purpose. Each prognostic equation can
be written in the general form

.. math::

   \label{eq:lfroggen}
       \frac{\partial\psi}{\partial t} + A^d(\psi)\psi + A^n(\psi)\psi = 0

where :math:`A^d` represents all adiabatic terms and :math:`A^n` all the
other terms (i.e., from parameterizations). Then the leapfrog scheme can
be given as

.. math::

   \label{eq:leapfrog}
       \frac{\psi^{t+\Delta t} - \psi^{t-\Delta t}}{2\Delta t} = -\left(A^d(\psi^t)\psi^t + A^n(\psi^{t-\Delta t})\psi^{t-\Delta t}\right)

where :math:`\Delta t > 0` is a predefined time increment. Note how in
equation `[eq:leapfrog] <#eq:leapfrog>`__ the adiabatic and nonadiabatic
terms are evaluated at different timesteps.

.. _prognostic-equations-1:

Prognostic Equations
~~~~~~~~~~~~~~~~~~~~

With these prerequisites at hand the discrete prognostic equations can
be given. The general procedure is straightforward: occurrences of
prognostic variables are replaced by their discrete values at the
respective grid points. Likewise, derivatives of these variables are
discretized according to the operator definitions given in the preceding
section. Wherever necessary average operators are used to ensure
consistency between all terms in the equation. All quantities in one
equation must either be located at mass points or at cell interfaces in
the horizontal and on either full or half-layers in the vertical.

Surface Pressure
^^^^^^^^^^^^^^^^

Approximating the integral in `[eq:ps] <#eq:ps>`__ with the midpoint
rule one has for the surface pressure at mass points

.. math::

   \label{eq:psdsc}
       \left(\frac{p_s^{t+\Delta t} - p_s^{t-\Delta t}}{2\Delta t}\right)_{ij} = -\frac{1}{a\cos{\varphi_j}}\sum_{l=1}^{k_m}{\left[\delta_\lambda U_l + \delta_\varphi\left(V_l\cos{\varphi}\right)\right]}

where :math:`k_m` is the number of vertical levels. Note that in REMO
:math:`k = 1` denotes the highest level and :math:`k = k_m` the lowest
one above the ground.

Wind Components
^^^^^^^^^^^^^^^

Equation `[eq:u] <#eq:u>`__ yields

.. math::

   \label{eq:udsc}
       \begin{split}
           &\left(\frac{u^{t+\Delta t}_k - u^{t-\Delta t}_k}{2\Delta t}\right)_{\iph j} - \frac{1}{\cos{\varphi_j}}\overline{Q_k}^\varphi\overline{V_k{\cos{\varphi}}}^{\lambda,\varphi} + \frac{1}{a\cos{\varphi_j}}\delta_\lambda\left(\Phi_k+K_k\right)\\ &+ \frac{R\overline{T_{vk}}^\lambda}{a\cos{\varphi_j}}\delta_\lambda\ln{p_k} + \frac{1}{\overline{\Delta p_k}^\lambda}\overline{\overline{\etastar}^\lambda\Delta_\eta\overline{u_k}^{2t}}^\eta = D_{uk} + F_{uk}
       \end{split}

for the horizontal wind component at :math:`u`-points. Likewise equation
`[eq:v] <#eq:v>`__ at :math:`v`-points results in:

.. math::

   \label{eq:vdsc}
       \begin{split}
           &\left(\frac{v^{t+\Delta t} - v^{t-\Delta t}}{2\Delta t}\right)_{i\jph} + \overline{Q_k}^\lambda\overline{U_k}^{\lambda,\varphi} + \frac{1}{a}\delta_\varphi\left(\Phi_k+K_k\right)\\ &+ \frac{R\overline{T_{vk}}^\varphi}{a}\delta_\varphi\ln{p_k}+ \frac{1}{\overline{\Delta p_k}^\varphi}\overline{\overline{\etastar}^\varphi\Delta_\eta\overline{v_k}^{2t}}^\eta = D_{vk} + F_{vk}
       \end{split}

Temperature
^^^^^^^^^^^

The discrete version of the thermodynamic equation `[eq:T] <#eq:T>`__
takes the form:

.. math::

   \label{eq:Tdsc}
       \begin{split}
           &\left(\frac{T_k^{t+\Delta t} - T_k^{t-\Delta t}}{2\Delta t}\right)_{ij} + \frac{1}{a\cos{\varphi_j}\Delta p_k}\left(\overline{U_k\delta_\lambda T_k}^\lambda + \overline{V_k\cos{\varphi}\delta_\varphi T_k}^\varphi\right)\\ &+ \frac{1}{\Delta p_k}\overline{\etastar\Delta_\eta\overline{T_k}^{2t}}^\eta = \alpha_k\omega_k + D_{Tk} + F_{Tk}
       \end{split}

Specific Humidity, Cloud Water and Cloud Ice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The structural similarity of equations
`[eq:qd] <#eq:qd>`__–`[eq:qi] <#eq:qi>`__ naturally carries over into
these discrete representations:

.. math::

   \label{eq:qddsc}
       \begin{split}
           &\left(\frac{q_{dk}^{t+\Delta t} - q_{dk}^{t-\Delta t}}{2\Delta t}\right)_{ij} + \frac{1}{a\cos{\varphi_j}\Delta p_k}\left(\overline{U_k\delta_\lambda q_{dk}}^\lambda + \overline{V_k\cos{\varphi}\delta_\varphi q_{dk}}^\varphi\right)\\ &+ \frac{1}{\Delta p_k}\overline{\etastar\Delta_\eta\overline{q_{dk}}^{2t}}^\eta = D_{q_dk} + F_{q_dk}
       \end{split}

.. math::

   \label{eq:qwdsc}
       \begin{split}
           &\left(\frac{q_{wk}^{t+\Delta t} - q_{wk}^{t-\Delta t}}{2\Delta t}\right)_{ij} + \frac{1}{a\cos{\varphi_j}\Delta p_k}\left(\overline{U_k\delta_\lambda q_{wk}}^\lambda + \overline{V_k\cos{\varphi}\delta_\varphi q_{wk}}^\varphi\right)\\ &+ \frac{1}{\Delta p_k}\overline{\etastar\Delta_\eta\overline{q_{wk}}^{2t}}^\eta = D_{q_wk} + F_{q_wk}
       \end{split}

.. math::

   \label{eq:qidsc}
       \begin{split}
           &\left(\frac{q_{ik}^{t+\Delta t} - q_{ik}^{t-\Delta t}}{2\Delta t}\right)_{ij} + \frac{1}{a\cos{\varphi_j}\Delta p_k}\left(\overline{U_k\delta_\lambda q_{ik}}^\lambda + \overline{V_k\cos{\varphi}\delta_\varphi q_{ik}}^\varphi\right)\\ &+ \frac{1}{\Delta p_k}\overline{\etastar\Delta_\eta\overline{q_{ik}}^{2t}}^\eta = D_{q_ik} + F_{q_ik}
       \end{split}

.. _diagnostic-equations-1:

Diagnostic Equations
~~~~~~~~~~~~~~~~~~~~

Equations `[eq:psdsc] <#eq:psdsc>`__–`[eq:qidsc] <#eq:qidsc>`__
constitute the discrete equation system that will be solved in each
timestep. The required discrete forms of the auxiliary and diagnostic
quantities will be given in the following.

.. _pressure-1:

Pressure
^^^^^^^^

Equation `[eq:akbkps] <#eq:akbkps>`__ yields the discrete pressure on
half layers:

.. math::

   \label{eq:pakbkdsc}
       p_{\kph} = A_{\kph} + B_{\kph}\cdot p_s

The discrete renditions of :math:`A` and :math:`B` here can simply be
derived by evaluating equations `[eq:ak] <#eq:ak>`__ and
`[eq:bk] <#eq:bk>`__ at the given discrete :math:`\eta`-levels. These
can as mentioned be arbitrarily chosen before the start of a simulation
and remain constant afterwards. Often pressure values on full layers
will be required instead:

.. math::

   \label{eq:pfdsc}
       p_k = \frac{p_\kph + p_\kmh}{2}

.. _potential-absolute-vorticity-1:

Potential Absolute Vorticity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equation `[eq:pav] <#eq:pav>`__ yields the discrete form

.. math::

   \label{eq:pavdsc}
       Q_k = \frac{1}{\overline{\Delta p_k a\cos{\varphi}}^{\lambda,\varphi}}\left[\overline{f}^{\lambda\varphi}a\overline{\cosphi}^\lambda + \delta_\lambda v_k - \delta_\varphi\left(u_k\cosphi\right)\right]

where discrete values of :math:`f` can be evaluated with equation
`[eq:coriolis] <#eq:coriolis>`__.

.. _geopotential-1:

Geopotential
^^^^^^^^^^^^

Applying once more the midpoint rule to the integral in equation
`[eq:geop] <#eq:geop>`__ gives for the geopotential on half layers:

.. math::

   \label{eq:geopdsch}
       \Phi_{\kph} = \Phi_s + R\sum_{l=k+1}^{k_m}{T_{vl}\ln{\frac{p_\lph}{p_\lmh}}}

From here full layer geopotentials are computed by averaging:

.. math::

   \label{eq:geopdscf}
       \Phi_k = \frac{\Phi_\kph + \Phi_\kmh}{2}

.. _kinetic-energy-1:

Kinetic Energy
^^^^^^^^^^^^^^

The discrete kinetic energy simply amounts to:

.. math::

   \label{eq:kindsc}
       K_k = \frac{1}{2}\left(\overline{u_k^2}^\lambda + \frac{1}{\cos{\varphi_j}}\overline{v_k^2\cosphi}^\varphi\right)

.. _virtual-temperature-1:

Virtual Temperature
^^^^^^^^^^^^^^^^^^^

For the discrete virtual temperature the discrete values of temperature
and moisture variables simply have to be substituted into equation
`[eq:vtemp] <#eq:vtemp>`__:

.. math::

   \label{eq:vtempdsc}
       T_{vk} = T_k\left(1 + \left(\frac{R_D}{R} - 1\right)q_{dk} - \left(q_{wk} + q_{ik}\right)\right)

.. _vertical-velocity-in-the-eta-system-1:

Vertical Velocity in the :math:`\eta`-System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Equation `[eq:etadot] <#eq:etadot>`__ reads in discrete form:

.. math::

   \label{eq:etadotdsc}
       \begin{split}
           &\etastar_{ij\kph} = b_\kph\frac{1}{a\cos{\varphi_j}}\sum_{l=1}^{k_m}{\left[\delta_\lambda U_l + \delta_\varphi\left(V_l\cos{\varphi}\right)\right]}\\ &-\frac{1}{a\cos{\varphi_j}}\sum_{l=1}^{k}{\left[\delta_\lambda U_l + \delta_\varphi\left(V_l\cos{\varphi}\right)\right]}
       \end{split}

Vertical Advection
^^^^^^^^^^^^^^^^^^

The vertical advection terms are discretized as follows

.. math::

   \begin{aligned}
       \label{eq:vadvudsc} 
       &\left(\etadot\dfdeta{u}\right)_{\iph jk} = \frac{1}{\overline{\Delta p_k}^\lambda}\overline{\overline{\etastar}^\lambda\Delta_\eta\overline{u_k}^{2t}}^\eta \\ &= \frac{1}{2\overline{\Delta p_k}^\lambda}\left[\overline{\etastar_\kph}^\lambda\left(\overline{u_{k+1}}^{2t} - \overline{u_k}^{2t}\right) + \overline{\etastar_\kmh}^\lambda\left(\overline{u_{k}}^{2t} - \overline{u_{k-1}}^{2t}\right)\right]\nonumber \\
       \label{eq:vadvvdsc}
       &\left(\etadot\dfdeta{v}\right)_{i\jph k} = \frac{1}{\overline{\Delta p_k}^\varphi}\overline{\overline{\etastar}^\varphi\Delta_\eta\overline{v_k}^{2t}}^\eta \\ &= \frac{1}{2\overline{\Delta p_k}^\varphi}\left[\overline{\etastar_\kph}^\varphi\left(\overline{v_{k+1}}^{2t} - \overline{v_k}^{2t}\right) + \overline{\etastar_\kmh}^\varphi\left(\overline{v_{k}}^{2t} - \overline{v_{k-1}}^{2t}\right)\right]\nonumber 
   \end{aligned}

where
:math:`\psi^{2t} = \frac{\psi^{t+\Delta t} - \psi^{t-\Delta t}}{2}`.
This approach is often referred to as the Crank-Nicolson scheme. The
important point is the presence of unknowns at time level
:math:`t + \Delta t` in the three different levels :math:`k-1, k` and
:math:`k+1`. This means that equations `[eq:vadvudsc] <#eq:vadvudsc>`__
and `[eq:vadvvdsc] <#eq:vadvvdsc>`__ constitute an implicit approach. As
will be seen later this requires the solution of a linear equation
system in every timestep.

.. _`subsec:solve`:

Solution of the Discrete Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the preceding sections the continuous equations
`[eq:ps] <#eq:ps>`__–`[eq:qi] <#eq:qi>`__ were rendered into the
discrete form `[eq:psdsc] <#eq:psdsc>`__–`[eq:qidsc] <#eq:qidsc>`__ to
facilitate the approximation of solutions on computer systems. In this
section it will be detailed how the resulting equation system can be
solved for the discrete values of the prognostic variables. With some
reordering all the discrete prognostic equations – with exception of the
surface pressure that does not contain implicit components – can be
written in the general form

.. math::

   \label{eq:abcd}
       A_k(\psi^{t+\Delta t}_{k-1}) + B_k(\psi^{t+\Delta t}_{k}) + C_k(\psi^{t+\Delta t}_{k+1}) = D_k

where :math:`\psi` is the prognostic variable of interest. For reasons
of brevity the exact form of the coefficients :math:`A_k, B_k, C_k` and
:math:`D_k` is not given here and the reader is instead referred to the
literature and specifically :raw-latex:`\textcite{emdm}`. Iterated over
all levels :math:`k=1,\dots,k_m` equation `[eq:abcd] <#eq:abcd>`__ then
yields a tridiagonal system that can be solved with standard methods. In
REMO the system is solved using the well-known Thomas algorithm. Note
that :math:`D_k` also contains values of the prognostic variable, albeit
only those at previous timesteps :math:`t-\Delta t` and :math:`t`.

Semi-Implicit Correction
~~~~~~~~~~~~~~~~~~~~~~~~

The previously presented scheme is – with the exception of vertical
advection – purely explicit. In practice however the model is almost
always run with additional implicit components. Since in REMO these
components can be switched on and off the term semi-implicit correction
is often used. First, the (mostly) explicit forecast is computed
according to the preceeding section. Then terms that are known to cause
gravity waves are corrected for by averaging the linear components of
those terms in time.

The main advantage of the semi-implicit approach is the ability to use
larger timesteps that would otherwise be prohibitive as a consequence of
the CFL-criterion :raw-latex:`\parencite[e.g.,][]{Duran2010}`. The
CFL-criterion essentially ensures that waves do not travel across more
than one gridbox per timestep. On the flipside the semi-implicit scheme
requires more computational effort, as the updated values of the
prognostic equations must first be derived from an appropriate equation
system. For REMO this results in a Helmholtz equation, that is quite
difficult to solve. Nevertheless the much larger timestep that can be
used with the semi-implicit correction generally outweighs this
disadvantage. This work will mostly be concerned with the explicit parts
of the model and hence a detailed description of the semi-implicit
correction is beyond the scope. Instead the reader is referred to the
literature and specifically to :raw-latex:`\cite{Simmons1981}` for an
in-depth treatment.

Asselin Filter
~~~~~~~~~~~~~~

To complete a prediction step the Asselin filter has to be applied to
the intermediate results. This is required to ensure the stability of
the leapfrog time integration scheme. With the definitions of equation
`[eq:lfroggen] <#eq:lfroggen>`__ the Asselin filter takes the form

.. math::

   \label{eq:asselin}
       \psi^{t+\Delta t} = \overline{\psi}^{t-\Delta t} - 2\Delta t A^d(\psi^t)\psi^t - 2\Delta t A^n(\psi^{t-\Delta t})\overline{\psi}^{t-\Delta t}

where :math:`\overline{\psi}` denotes the filtered values.

This completes the introduction into the REMO dynamical core. In the
next chapter the notations and definitions given here will serve as the
basis for taking a closer look at the pressure gradient force in REMO
and the particular problems arising from its computation.

.. |image1| image:: images/rotation.png
.. |image2| image:: images/discretization.png
.. |image3| image:: images/cgrid.png
.. |image4| image:: images/lorenz.png
