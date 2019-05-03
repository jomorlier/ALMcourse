# From Finite Element Analysis to Topology Optimization

This note introduces the basic concept of finite element analysis (FEA) and topology optimization (Topopt), 
and hopefully convince you that these are important topics to learn about for future 
engineers 

We showed that the deflection or deformation of a structure (${\bf U}$) can be solved through a algebraic equation 

$${\bf KU}={\bf F}$$.

It can be further shown that the strain energy due to the deformation is $$0.5{\bf U}^T{\bf KU}$$. Minimizing this energy
 is equivalent to minimizing the compliance of the structure under the given loads and boundary conditions, leading to 
 optimal topologies. 

We explain some technical details below. A good tutorial can be found from Dr. Sigmund's group (see e.g., 
this [code](http://www.topopt.dtu.dk/?q=node/751) and this [paper](http://www.topopt.dtu.dk/files/TopOpt88.pdf)).

### The compliance minimization problem

Topology optimization has been commonly used to design structures and materials with optimal mechanical, 
thermal, electromagnetic, acoustical, or other properties. 
The structure under design is segmented into $$n$$ finite elements, and a 
density value $$x_i$$ is assigned to each element $$i \in \{1,2,...,n\}$$: A higher density corresponds to a less porous material 
element and higher Yong's modulus. Reducing the density to zero is equivalent to creating a hole in the 
structure. Thus, the set of densities $${\bf x}=\{x_i\}$$ can be used to represent the topology of the 
structure and is considered as the variables to be optimized. A common topology optimization problem is 
compliance minimization, where we seek the "stiffest" structure within a certain volume limit to withhold 
a particular load: 

$$ \min_{\bf x} \quad {\bf f} := {\bf U}^T {\bf K}({\bf x}) {\bf U} $$

$$ \text{subject to:} \quad {\bf h} := {\bf K}({\bf x}) {\bf U} = {\bf F}, $$

$$ \quad {\bf g} := V(\textbf{x}) \leq v,$$

$$ \textbf{x} \in [0,1]. $$

Here $$V(\textbf{x})$$ is the total volume; $$v$$ is an upper bound on volume; 
$${\bf U} \in \mathbb{R}^{n_d\times 1}$$ is the displacement of the structure under the load $${\bf F}$$, 
where $$n_d$$ is the degrees of freedom (DOF) of the system (i.e., the number of x- and y-coordinates 
of nodes from the finite element model of the structure);
$${\bf K(x)}$$ is the global stiffness matrix for the structure.

$${\bf K(x)}$$ is indirectly influenced by the topology $${\bf x}$$, through the element-wise stiffness matrix 

$${\bf K}_i = \bar{\bf K}_e E(x_i)$$,

where the matrix $$\bar{\bf K}_e$$ is predefined according to the finite element type 
and the nodal displacements of the element, $$E(x_i)$$ is 
the element-wise Young's modulus defined as a function of the density $$x_i$$: 
$$E(x_i) := \Delta E x_i^p + E_{\text{min}}$$, where $$p$$ (the penalty parameter) is usually set to 3.
 This cubic relationship between the topology and the modulus 
is determined by the material constitutive models, and numerically, it 
also helps binarize the topologies, i.e., to push the optimal $${\bf x}_i$$ to 1 or 0 (why?). 
The term $$E_{\text{min}}$$ is added to provide numerical stability.

For those interested, [here](http://designinformaticslab.github.io/designopt_tutorial/2017/10/26/topologyopt.html) 
I explain the details for solving this topology optimization problem.

## Application of Topology Optimization

While the theories behind topology optimization were developed during the 1980s, the value of this technique 
only started to show after manufacturing and in particular 3D printing technologies become 
more matured in the last two decades. Today topology optimization has many industry applications, from 
light-weight vehicle bodies, to artificial bones and implants, to biomemetic airplane wings, 
to various sensors and actuators.

<img src="/_images/mechdesign/featop4.jpg" alt="Drawing" style="height: 300px;"/> 

<img src="/_images/mechdesign/featop7.png" alt="Drawing" style="height: 300px;"/> 

Image from Aage et al. (2017) "Giga-voxel computational morphogenesis for structural design"

<img src="/_images/mechdesign/featop5.png" alt="Drawing" style="height: 300px;"/> 

Image from Wu et al. (2018) "Infill Optimization for Additive Manufacturing - Approaching Bone-like Porous Structures"

<img src="/_images/mechdesign/featop6.png" alt="Drawing" style="height: 400px;"/> 

Image from Dassault Systems

You can find in teaching repository 

[The course introduction](https://github.com/jomorlier/ALMcourse/blob/master/teaching/TOPOPT_intro2018.pdf)

[The exercice using top88](https://github.com/jomorlier/ALMcourse/blob/master/teaching/BE_Topopt_eleve.pdf)

[The highlights for FA & TOPOPT](https://github.com/jomorlier/ALMcourse/blob/master/teaching/Higlights_FA_JM.pdf)

[Our TOPOPT+ALM Research](https://github.com/jomorlier/ALMcourse/blob/master/teaching/TOPOPT&ALM_SUPAERO.pdf)


[The 3 point bending projected corrected using top88](http://htmlpreview.github.io/?https://github.com/jomorlier/ALMcourse/blob/master/top88/topopt_3ptBENDING.html)


