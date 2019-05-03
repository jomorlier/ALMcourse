# crash course

You can find in teaching repository 

[The course introduction](https://github.com/jomorlier/ALMcourse/blob/master/teaching/TOPOPT_intro2018.pdf)

[The exercice using top88](https://github.com/jomorlier/ALMcourse/blob/master/teaching/BE_Topopt_eleve.pdf)

[The highlights for FA & TOPOPT](https://github.com/jomorlier/ALMcourse/blob/master/teaching/Higlights_FA_JM.pdf)

[Our TOPOPT+ALM Research](https://github.com/jomorlier/ALMcourse/blob/master/teaching/TOPOPT&ALM_SUPAERO.pdf)

[The 3 point bending projected corrected using top88](https://github.com/jomorlier/ALMcourse/blob/master/top88/topopt_3ptBENDING.html)



# From Finite Element Analysis to Topology Optimization

This note introduces the basic concept of finite element analysis (FEA) and topology optimization (Topopt), 
and hopefully convince you that these are important topics to learn about for future 
engineers 

We showed that the deflection or deformation of a structure ($${\bf U}$$) can be solved through a algebraic equation 

$${\bf KU}={\bf F}$$.

It can be further shown that the strain energy due to the deformation is $$0.5{\bf U}^T{\bf KU}$$. Minimizing this energy
 is equivalent to minimizing the compliance of the structure under the given loads and boundary conditions, leading to 
 optimal topologies. 

We explain some technical details below. A good tutorial can be found from Dr. Sigmund's group (see e.g., 
this [code](http://www.topopt.dtu.dk/?q=node/751) and this [paper](http://www.topopt.dtu.dk/files/TopOpt88.pdf)).


For those interested, [here](http://designinformaticslab.github.io/designopt_tutorial/2017/10/26/topologyopt.html) 
I explain the details for solving this topology optimization problem.

## Application of Topology Optimization

While the theories behind topology optimization were developed during the 1980s, the value of this technique 
only started to show after manufacturing and in particular 3D printing technologies become 
more matured in the last two decades. Today topology optimization has many industry applications, from 
light-weight vehicle bodies, to artificial bones and implants, to biomemetic airplane wings, 
to various sensors and actuators.

<img src="/_images/mechdesign/featop4.jpg" alt="Drawing" style="height: 150px;"/> 

<img src="/_images/mechdesign/featop7.png" alt="Drawing" style="height: 150px;"/> 

Image from Aage et al. (2017) "Giga-voxel computational morphogenesis for structural design"

<img src="/_images/mechdesign/featop5.png" alt="Drawing" style="height: 1500px;"/> 

Image from Wu et al. (2018) "Infill Optimization for Additive Manufacturing - Approaching Bone-like Porous Structures"

<img src="/_images/mechdesign/featop6.png" alt="Drawing" style="height: 200px;"/> 

Image from Dassault Systems



