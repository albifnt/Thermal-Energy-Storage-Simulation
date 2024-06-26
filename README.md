# Thermal-Energy-Storage-Simulation

## Simulation of a Thermal Energy Storage using Computational Fluid Dynamics methods
<p align="justify">
The availability of energy has always represented a fundamental need for the
human life, and it will continue to be like this in the future as well. In the past energy
was mainly provided by fossil fuels, but because of their shortage and the interest in
cleaner and more efficient energy sources, engineers started to look for new ways of
producing energy. The energy coming from the sun is the most plentiful source human
being can exploit on the Earth, with 885 milion of Terawatt hours (TWh) reaching the
terrestrial surface. The two technologies which are mainly used take advantage of solar
energy are the Photovoltaic System and the Concentrated Solar Power (CSP).
A CSP system uses the solar radiation as energy source at high temperature. In
order to reach temperatures which are competitive with those reached through fossil
fuels, the radiation is concentrated on a solar receptacle, where it is then transferred to
the heat transfer fluid as thermal energy. The fluid flows along the receptacles and
transfers the heat, which, in case of production of electric energy, it is converted into
superheated steam that activates the turbine through a thermodynamic cycle.
</p>

<p align="center">
  <img  style="text-align:center" src="https://github.com/albifnt/Thermal-Energy-Storage-Simulation/blob/main/images/Figure_1.JPG">
  <br>
  Fig.1 - Basic operation of a CSP.
</p>

<p align="justify">
As shown in Fig.1, CSPs cannot perform their activity when solar radiation does not
hit the Earth surface. Anyway, these systems present the intrinsic characteristic of
storing energy, thanks to their thermal inertia, and this is what makes this technology
really interesting. In the CSP technology, the word “thermal storage” indicates the
accumulation of the exceeding thermal energy during the daily hours, so that it can be
available after the sunset or in other moments of the daytime so as to supply all the
needs required by the electric network.
The exceeding energy is stored in the so called Thermal Energy Storage (TES). This
system is active if the storage medium can move, while it is referred to as passive if the
storage medium is stationary. The code in this repository simulates a passive TES. In such a
system, hot and cold temperatures coexist, yielding to a vertical temperature gradient.
The heat transfer fluid transports energy to the medium during the charging phase, and
receives it from the same medium during the discharging phase, taking it to the turbine
in those moments of the day when the solar radiation is missing.
</p>

<p align="center">
  <img src="https://github.com/albifnt/Thermal-Energy-Storage-Simulation/blob/main/images/Figure_2.JPG">
  <br>
  Fig.2 - Basic operation of a TES.
</p>

## Governing equations
<p align="justify">
The equations governing the temparature exchange between the fluid and the solid are:
</p>

<p align="center">
$` ε ρ_f C_{p,f} {∂T_f\over ∂t} + ε ρ_f C_{p,f} u_{f,i} {∂T_f\over ∂x} = k_f {∂^2T_f\over ∂x^2} + h_f ( T_s - T_f ) `$
</p>
<p align="center">
$` (1 - ε) ρ_s C_s {∂T_s\over ∂t} = k_s {∂^2T_s\over ∂x^2} - h_v ( T_s - T_f ) `$
</p>

where:

$` ε `$ : porosity

$` ρ_f `$ : fluid density

$` ρ_s `$ : solid density

$` C_{p,f} `$ : fluid specific heat at constant pressure

$` C_s `$ : solid specific heat

$` k_f `$ : fluid thermal conductivity

$` k_s `$ : solid thermal conductivity

$` h_f `$ : fluid heat transfer

$` h_v `$ : volumetric heat-transfer

<p align="justify">
The numerical method adopted for the discretization of the governing equation is
the Finite Volume Method (FVM). The FVM is based on two basic ideas:
  
1. Solving the integral form of the governing equation
   
2. Evaluating fluxes through the cell averages by using “reconstruction”
</p>

## Order Verification Study
<p align="justify">
The purpose of the Order Verification Study consists of analysing and evaluating
the spatial discretization of the approximate equations (only the spatial discretization).
This analysis is extremely important, precise and it also offers the possibility of finding
programming mistakes. By testing the implementation through the Order Verification
Study, we can determine if the spatially discretized equations are consistent or not: a
model achieves consistency when it approaches to the exact solution as the grid spacing
is reduced, while it is called inconsistent when this does not happen.
According to the definition of consistency, the implementation of the equation is tested
by reducing the grid spacing and monitoring if the results approaches to the exact
solution. This is called Grid-Refinement Study. The only problem is that we do not
know the exact solution of the mathematical model, for this reason we need a method that allows
to assess the consistency of the semi-discretized model and at the same time to check if
its actual order of accuracy is equal to the nominal order of accuracy. This can be achieved by 
using the Method of Manufactured Solution (MMS). The idea behind the MMS is really simple and can be
explained in 5 points:
  
1. Choice of the form of the manufactured solution
   
2. Generation of the analytical source term that is obtained by applying the
governing equation to the manufactured solution

3. Resolution of the discretized governing equation with the analytical
result got at point 2 acting as source term, by reducing the grid spacing

4. Evaluation of the errors between the exact manufactured solution and the
discretized manufactured solution

5. Determination of the order of accuracy
</p>

## Storage Design Study
<p align="justify">
The Storage Design Study consists of finding the optimal dimensions (diameter and height) given 
the volume, in other words, the dimensions which maximise the exergy efficiency and the capacity
factor. The exergy efficiency is defined as:
</p>

<p align="center">
$` \frac{Ex_{d,out} - Ex_{d,in}}{Ex_{c,in} - Ex_{c,out}} `$
</p>

<p align="justify">
where the symbol ̇<b>Ex</b> is the exergy flux (the subscripts <b>c, d</b> indicate the charging and the
discharging phase respectively and <b>in, out</b> indicate the inlet and the outlet), defined as:
</p>

<p align="center">
$` ∫ mx_f C_{p,f} [T - T_0 - T_0   ln(\frac{T_0}{T})] dt  `$
</p>

<p align="justify">
where <b>mx_f</b> is the mass flow rate.
</p>

## Structure of the code
<p align="justify">
The code is composed of three main parts:
  
1. The Initialization

2. The Main

3. The Definition of the functions


The first part is the one dedicated to the initialization of functions, variables and
files which contain the data necessary for plots and graphs. In the main, the user can
take advantage of the different functions to determine the status of the system, to carry
on the Order Verification Study and to conduct the Storage Design Study. Let’s now
describe how to generate the different simulations:
</p>

#### The Order Verification Study
<p align="justify">
In order to carry on the Order Verification Study it is necessary to call the
function <b>Order_Verification_Study()</b> inside the main of the code, by
commenting the other functions. If the user wants to run the OVS for different
values of the paramaters, he should change these last ones in the function
<b>input_OVS()</b>, leaving the variables <b>flag=1</b> and <b>DOS=0</b> unchanged, since they
are essential to implement just the OVS. If in the OVS, the user wants to
decouple the fluid and solid equations, he should impose the variables <b>hv_f</b> and
<b>hv_s</b> equal to zero in <b>input_OVS()</b>.
</p>

#### Exact solution
<p align="justify">
To carry on the comparison with the exact solution, it is
necessary to call the function <b>Simulation()</b> inside the main of the code, by
commenting the other functions. If the user wants to run <b>Simulation()</b> for
different values of the paramaters, he should change these last ones in the
function <b>input_Simulation()</b>, leaving the variables <b>flag=0</b> and <b>DOS=0</b>
unchanged, since they are necessary to implement just the Exact solution.
</p>

#### Storage Design Study
<p align="justify">
To carry on the Storage Design Study, it is necessary to
call the function <b>Storage_Design_Study()</b> inside the main of the code, by
commenting the other functions. If the user wants to run
<b>Storage_Design_Study()</b> for different values of the paramaters, he should
change these last ones in the function <b>input_Storage_Design_Study()</b>,
leaving the variables <b>flag=0</b> and <b>DOS=1</b> unchanged, since they are necessary
to implement just the Storage Design Study.
</p>

## MATLAB scripts
<p align="justify">
The graphs can be plotted using the different Matlab scripts:
</p>

#### Exact_solution
<p align="justify">
This Matlab script can be used to draw the thermocline of the charging phase in order to compare it with the exact solution. The graph contains the plot of the exact solution and another one obtained from a simulation with n_cells=50. Before running this script, it is necessary that the user has run <b>Simulation()</b>.
</p>

#### Fluid_GRS – Solid_GRS
<p align="justify">
This two Matlab scripts are used to carry on the Grid Refinement Study for the solid and the fluid. They plot the graphs of the Error and of the actual order of accuracy as functions of the grid spacing. Furthermore, a graph shows how the manufactured solution approximate the exact solution. Before running these scripts, it is necessary that the user has run <b>Order_Verification_Study()</b>.
</p>

#### Storage_Design_Study
<p align="justify">
This Matlab script allows to draw Exergy Efficiency, the Capacity Factor and the Temperature increase from the Storage Design Study. Before running this script, it is necessary that the user has run <b>Storage_Design_Study</b>.
</p>

#### Temperature_evolution
<p align="justify">
This Matlab script plots a moving graph of the fluid and solid temperature throughout the different phases of a whole cycle.
</p>
