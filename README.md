<font color="blue"><b><h1 style="font-size:4em; font-family:serif"><center>GALACTIC DYNAMO</center></h1></b>

This project is a part of the term project given in course P464 Plamsa Physics and Magnetohydrodynamics taught in Spring 2024 at NISER Bhubaneswar.

Submitted by: <u>Chandan Kumar Sahu</u>, Integrated MSc. SPS batch 19

Supervised by: Dr. Luke R. Chamandy, SPS, NISER

# Contents
1. Problem Statement
2. Theory
    1. The mean-field induction equation
    2. Solving in cylindrical coordinates
    3. Calculation of total magnetic field magnitude and pitch angle
    4. Magnetic decay constant
4. Numerical solution to Galactic Magnetic Fields
    1. Solver: Crank-Nicholson
    2. Defining the grid 
    3. Boundary conditions

    # 1 - Problem Statement

Solve the diffusion equation in the z-direction.
* Explore the evolution of the magnetic field magnitude and of the exponential decay rate.
* Explore the evolution of the spatial solution for $B_r$ and $B_\phi$, and of the pitch angle of the mean magnetic field p.
* Explore how different boundary conditions affect the results.
* Explore how different seed fields affect the results.

<font color="blue"><b><h1 style="font-size:4em; font-family:serif"><center>GALACTIC DYNAMO</center></h1></b>

This project is a part of the term project given in course P464 Plamsa Physics and Magnetohydrodynamics taught in Spring 2024 at NISER Bhubaneswar.

Submitted by: <u>Chandan Kumar Sahu</u>, Integrated MSc. SPS batch 19

Supervised by: Dr. Luke R. Chamandy, SPS, NISER


# Contents
1. Problem Statement
2. Theory
    1. The mean-field induction equation
    2. Solving in cylindrical coordinates
    3. Calculation of total magnetic field magnitude and pitch angle
    4. Magnetic decay constant
4. Numerical solution to Galactic Magnetic Fields
    1. Solver: Crank-Nicholson
    2. Defining the grid 
    3. Boundary conditions
# 1 - Problem Statement

Solve the diffusion equation in the z-direction.
* Explore the evolution of the magnetic field magnitude and of the exponential decay rate.
* Explore the evolution of the spatial solution for $B_r$ and $B_\phi$, and of the pitch angle of the mean magnetic field p.
* Explore how different boundary conditions affect the results.
* Explore how different seed fields affect the results.
# 2 - Theory

The dynamo theory describes the process through which a rotating, convecting, and electrically conducting fluid can maintain a magnetic field over astronomical time scales.
The evolution of the magnetic field is described by the induction equation, which relates changes in the magnetic field to the velocity field of the conducting fluid and its magnetic diffusivity. 

In this article, we will describe how the galactic magnetic field varies with time in the $z$-direction, and discuss the points mentioned in the probelm statement.

## 2.1 - The mean-field induction equation

We have the mean-field induction equation as 
$$ \dfrac{\partial \bar{\mathbf{B}}}{\partial t} = \nabla \times \left[ \bar{\mathbf{V}} \times \bar{\mathbf{B}} + \mathcal{E} - \eta \left( \nabla \times \bar{\mathbf{B}} \right) \right] $$
where $\mathcal{E} = \left( \alpha \bar{\mathbf{B}} \right) - \eta_t \left( \nabla \times \bar{\mathbf{B}} \right)$

We will solve the equations in the cylindrical coordinates (r, $\phi$, z) with the origin at the galactic centre and the z-axis parallel to the galactic angular velocity. However, to simplify things, lets make some approximations.

1. Omit the terms involving $\bar{\mathbf{V}} \times \bar{\mathbf{B}}$ and $\alpha$. We will land up on just the diffusion equation.
$$ \dfrac{\partial \bar{\mathbf{B}}}{\partial t} = - \nabla \times \left[ \eta_T \left( \nabla \times \bar{\mathbf{B}} \right) \right] $$
where $\eta_T = \eta + \eta_t$

2. Take $\eta_T$ independent of $\bar{\mathbf{B}}$. Our equation becomes
$$ \dfrac{\partial \bar{\mathbf{B}}}{\partial t} = - \eta_T \left[ \nabla \times \left( \nabla \times \bar{\mathbf{B}} \right) \right] $$
But $\nabla \times \left( \nabla \times \bar{\mathbf{B}} \right) = \nabla \left( \nabla \cdot \bar{\mathbf{B}} \right) - \nabla^2 \bar{\mathbf{B}} $ and $\nabla \cdot \bar{\mathbf{B}} = 0$ (Gauss's Law), so we finally have
$$ \boxed{ \dfrac{\partial \bar{\mathbf{B}}}{\partial t} = \eta_T \nabla^2 \bar{\mathbf{B}} }$$

This is the Fickian diffusion equation. We will solve this equation numerically.

## 2.2 - Solving in cylindrical coordinates

In cylindrical coordinates, 
$$ \mathbf{\bar{B}} = \bar{B}_r \mathbf{\hat{r}} + \bar{B}_{\phi} \mathbf{\hat{\phi}} + \bar{B}_z \mathbf{\hat{z}} $$

The $\nabla^2B$ becomes
$$ \begin{aligned}
\nabla^2 \mathbf{\bar{B}} = & \left[ \frac{\partial}{\partial r} \left( \frac{1}{r} \frac{\partial}{\partial r} \left( r \bar{B}_r \right) \right) + \frac{1}{r^2} \frac{\partial^2 \bar{B}_r}{\partial \phi^2} + \frac{\partial^2 \bar{B}_r}{\partial z^2} - \frac{2}{r^2} \frac{\partial \bar{B}_\phi}{\partial \phi} \right] \mathbf{\hat{r}} \\
& + \left[ \frac{\partial}{\partial r} \left( \frac{1}{r} \frac{\partial}{\partial r} \left( r \bar{B}_\phi \right) \right) + \frac{1}{r^2} \frac{\partial^2 \bar{B}_\phi}{\partial \phi^2} + \frac{\partial^2 \bar{B}_\phi}{\partial z^2}+\frac{2}{r^2} \frac{\partial \bar{B}_r}{\partial \phi} \right] \mathbf{\hat{\phi}} \\
& + \left[\frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial \bar{B}_z}{\partial r}\right)+\frac{1}{r^2} \frac{\partial^2 \bar{B}_z}{\partial \phi^2}+\frac{\partial^2 \bar{B}_z}{\partial z^2}\right] \mathbf{\hat{z}} 
\end{aligned} $$

So we get the final component-wise equations as
$$ \begin{aligned}
    \frac{\partial \bar{B}_r}{\partial t} &= \eta_T \left[ \frac{\partial}{\partial r} \left( \frac{1}{r} \frac{\partial}{\partial r} \left( r \bar{B}_r \right) \right) + \frac{1}{r^2} \frac{\partial^2 \bar{B}_r}{\partial \phi^2} + \frac{\partial^2 \bar{B}_r}{\partial z^2} - \frac{2}{r^2} \frac{\partial \bar{B}_\phi}{\partial \phi} \right] \\
    \frac{\partial \bar{B}_\phi}{\partial t} &= \eta_T \left[ \frac{\partial}{\partial r} \left( \frac{1}{r} \frac{\partial}{\partial r} \left( r \bar{B}_\phi \right) \right) + \frac{1}{r^2} \frac{\partial^2 \bar{B}_\phi}{\partial \phi^2} + \frac{\partial^2 \bar{B}_\phi}{\partial z^2}+\frac{2}{r^2} \frac{\partial \bar{B}_r}{\partial \phi} \right] \\
    \frac{\partial \bar{B}_z}{\partial t} &= \eta_T \left[\frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial \bar{B}_z}{\partial r}\right)+\frac{1}{r^2} \frac{\partial^2 \bar{B}_z}{\partial \phi^2}+\frac{\partial^2 \bar{B}_z}{\partial z^2}\right] 
\end{aligned} $$

Since the problem statement in the project aims to solve only in the z-direction, we remove all radial or azimuthal variations of the magnetic field $\left(\dfrac{\partial }{\partial r} = \dfrac{\partial }{\partial \phi} = 0 \right)$. 

We are now left with these simple equations to solve, i.e., the Fickian diffusion equations.
$$ \boxed{ \frac{\partial \bar{B}_r}{\partial t} = \eta_T \frac{\partial^2 \bar{B}_r}{\partial z^2} } \qquad \qquad \qquad \boxed{ \frac{\partial \bar{B}_\phi}{\partial t} = \eta_T \frac{\partial^2 \bar{B}_\phi}{\partial z^2} } \qquad \qquad \qquad \boxed{ \frac{\partial \bar{B}_z}{\partial t} = \eta_T \frac{\partial^2 \bar{B}_z}{\partial z^2} } $$


## 2.3 - Calculation of total magnetic field magnitude ($B_{total}$) and pitch angle ($p_B$)

We can calculate the magnitude of the total magnetic field as 
$$ B_{\text{total}} = \sqrt{\bar{B}_r^2 + \bar{B}_\phi^2} $$
And $$ p_B = \tan^{-1} \left( \dfrac{\bar{B}_r}{\bar{B}_\phi} \right) \qquad \qquad \text{where} \qquad -\dfrac{\pi}{2} < p_B < \dfrac{\pi}{2}$$
## 2.4 - Magnetic decay constant

The diffusion of the total magnetic field can be expressed in the form

$$ B_{\text{total}}(z, t) = \tilde{B}(r)\exp(-\gamma t) $$ 
where $\tilde{B}(r)$ contains all the variation in $r$ and the exponential factor contains time variation, where $\gamma$ is the magnetic decay constant.

Our final goal is to calculate this $\gamma$.
# 3 - Numerical Solution to Galactic magnetic Fields


Now we begin to solve the diffusion equations numerically.

We will first solve the heat diffusion equation, and use them to calculate the total $B_{\text{total}}$ and the pitch angle $p_B$. Finally we will calculate the decay constant $\gamma$.

## 3.1 - Solver: Crank-Nicholson
We use the Crank-Nicholson algorithm to solve the diffusion equation. It combines implicit and explicit schemes, resulting in a stable and accurate solution. 

A stencil of forward eluer method and crank-nicholson method are shown below:

Forward Euler Method

![euler.png](attachment:euler.png)


Crank Nicholson Method

![cn.png](attachment:cn.png)

By averaging the values of variables at current and next time steps, it reduces numerical errors and suppresses oscillations. 

The method is highly stable, allowing for larger time steps compared to explicit methods, and is second-order accurate in time and space. This makes it well-suited for simulating heat diffusion processes with high accuracy and efficiency.

## 3.2 - Defining the grid

#### $\star$  Spatial Grid ($z$)

The size of a galactic disk typically ranges from around 10 kpc to 50 kpc. Our Milky way has a radius of ~40 kpc. Within this disk, there are thin and thick disks. The thin disk, where most of the younger stars reside, typically extends from the galactic center to about 100-1000 pc. Meanwhile, the thick disk, containing older stars and having a higher vertical velocity dispersion, can extend further out to ~3 kpc. 

Considering these typical values, we choose the thickness of our model galaxy to be ~200 pc. 

So, $z$ extends from -100 pc to +100 pc. For solving numerically, we normalize the spatial grid to -1 to 1, with a grid cell size of $dz = 0.01$.

![Milky-way.jpg](attachment:Milky-way.jpg)
<br>
<font>Figure: Schematic of the Milky Way galaxy

#### $\star$  Temporal Grid ($t$)

The typical variation of magnetic fields usually vary in millions of years. 
However, the duration of simulation depends on the seed field, so we change this values and its binning according the the seed field. In general, we take 100 points within the temporal grid range.

#### $\star$  Magnetic diffusion coefficient ($\eta_T$)

The magnetic diffusion coefficient measures the rate at which magnetic fields diffuse through a medium. We calculate the value of $\eta_T$ as
$$ \eta_T \approx \dfrac{1}{3} \tau v_{\text{rms}}^2 $$

But we first need to calculate its value in the normlizated units.
For a typical galaxy, the faraday time $\tau \approx 10$ Myr and velocity $v_{\text{rms}} \approx 10^4$ km/s

$$ \begin{aligned}
\eta_T &= \dfrac{1}{3} \times 10 \text{ Myrs} \times (10 \text{ km/s})^2 \\
&= \dfrac{100}{3} \times \text{ Myrs} \times  \dfrac{\text{ km}^2}{\text{s}^2} \\
&= \dfrac{100}{3} \times \text{ Myrs} \times  \left( \dfrac{\text{ km} \times 100 \text{ pc}}{3.086\times 10^{15}\text{ km}} \right)^2 \left(\dfrac{3.154\times 10^{13} \text{ s}}{\text{s} \times \text{Myr}} \right)^2 \\
&= \dfrac{100}{3} \times \left( \dfrac{3.154\times 10^{13}}{3.086\times 10^{15}} \right)^2 \times \text{ (100 pc)}^2/\text{ Myr} \\
&= 3.48 \times 10^{-3}\text{ (100 pc)}^2/\text{ Myr}
\end{aligned} $$

Hence, we use the following values
$$ -1<z<1 \qquad \qquad \text{ with } dz=0.01 \text{ and } \qquad \qquad \eta_T = 3.48 \times 10^{-3}\text{ (100 pc)}^2/\text{ Myr} $$

## 3.3 - Boundary conditions

Two common types of boundary conditions are vacuum boundary conditions and isolated boundary conditions.

#### Vacuum Boundary Conditions (BC):
In vacuum boundary conditions, the system is assumed to be in contact with a perfect vacuum, the boundaries of the system under consideration are devoid of any matter or energy. It ensures that the quantity becomes 0 at the boundaries at all time.

#### Isolated Boundary Conditions:
Isolated boundary conditions imply that the system is completely isolated from its surroundings, meaning there are no exchanges of matter or energy across its boundaries. It ensures that the derivative of the quantity becomes 0 at the boundaries at all time.

We are calculating the diffusion of $B_r$ and $B_\phi$. As the distance in z-increases beyond, $B_z$ dominates over $B_r$. 

Within the thin disk approximation, it is safe to assume that $B_r \rightarrow 0$ as $z \rightarrow 100$ pc. The same thing applies to $B_\phi$ too. Hence, we apply vacuum boundary conditions for all the cases here.
