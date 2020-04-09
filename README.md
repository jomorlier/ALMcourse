# TOPOPT 2020

****

LECTURES

****

[The course introduction](https://github.com/jomorlier/ALMcourse/blob/master/teaching/TOPOPT_intro2018.pdf)

with audio using Flash

[Flash](https://www.commentcamarche.net/contents/447-swf-fichier-swf)

[Audio1](https://github.com/jomorlier/ALMcourse/blob/master/TOpopt1.swf)

[Audio2](https://github.com/jomorlier/ALMcourse/blob/master/TOpopt2.swf)

[Audio3](https://github.com/jomorlier/ALMcourse/blob/master/TOpopt3.swf)


[But how do you compute sensitivity of the compliance?](https://github.com/jomorlier/ALMcourse/blob/master/teaching/C3_demo.pdf)

****

Jack Chessa's Youtube channel

[Hypermesh/hyperview tutorials](https://www.youtube.com/playlist?list=PL3A7B78F0E428DF72)

****



****

COMPUTER LAB using SIMP

****

[The exercice using top88](https://github.com/jomorlier/ALMcourse/blob/master/teaching/BE_Topopt_eleve.pdf)

[The highlights for FA & TOPOPT](https://github.com/jomorlier/ALMcourse/blob/master/teaching/Higlights_FA_JM.pdf)

[Our TOPOPT+ALM Research](https://github.com/jomorlier/ALMcourse/blob/master/teaching/TOPOPT&ALM_SUPAERO.pdf)

You can find in top88 repository all the files, paper, and ...

[The 3 point bending projected corrected using top88](http://htmlpreview.github.io/?https://github.com/jomorlier/ALMcourse/blob/master/top88/topopt_3ptBENDING.html)

for people wondering about how linear elasticity of top88.m is working

Prof, can you help to us understand how to compute the stiffness matrix of 2D membrane?
Explictely ?  [Membrane2D_K](http://htmlpreview.github.io/?https://github.com/jomorlier/feacourse/blob/master/Membrane2D_K/Elementarystiffrecmesh.html)


Have a look in the top3D repository to do the same but in 3D !!! 


****

MORE ADVANCED STUFF

****

Last but not the least a MATLAB tutorial Advanced Topology Optimization.

Part A:  [Constraints Agreggation](http://htmlpreview.github.io/?https://github.com/jomorlier/ALMcourse/blob/master/AdvancedTopOpt/ConstraintsAgreggation.html)
Thanks to my PhD Simone.

Part B:  [Stress Based TopOpt](http://htmlpreview.github.io/?https://github.com/jomorlier/ALMcourse/blob/master/AdvancedTopOpt/StressBasedTopOpt.html)
Thanks **AGAIN** to my PhD Simone.
Before you can use Method of Moving Asymptotes (MMA) as an optimizer in our stress based topology optimization program, you need to obtain the Matlab implementation of MMA from Prof. Krister Svanberg (krille@math.kth.se) from KTH in Stockholm Sweden.



****

RECAP

****
# From Finite Element Analysis to Topology Optimization

This note introduces the basic concept of finite element analysis (FEA) and topology optimization (Topopt), 
and hopefully convince you that these are important topics to learn about for future 
engineers 

We showed that the deflection or deformation of a structure U can be solved through a algebraic equation 

KU=F.

It can be further shown that the strain energy due to the deformation is 0.5U^TKU. Minimizing this energy
 is equivalent to minimizing the compliance of the structure under the given loads and boundary conditions, leading to 
 optimal topologies. 

We explain some technical details below. A good tutorial can be found from Dr. Sigmund's group (see e.g., 
this [code](http://www.topopt.dtu.dk/?q=node/751) and this [paper](http://www.topopt.dtu.dk/files/TopOpt88.pdf)).


For those interested, A good introduction can be find [here](http://designinformaticslab.github.io/designopt_tutorial/2017/10/26/topologyopt.html) 


## Application of Topology Optimization

While the theories behind topology optimization were developed during the 1980s, the value of this technique 
only started to show after manufacturing and in particular 3D printing technologies become 
more matured in the last two decades. Today topology optimization has many industry applications, from 
light-weight vehicle bodies, to artificial bones and implants, to biomemetic airplane wings, 
to various sensors and actuators.

<img src="/_images/mechdesign/featop4.jpg" alt="Drawing" style="height: 100px;"/> 

<img src="/_images/mechdesign/featop7.png" alt="Drawing" style="height: 100px;"/> 

Image from Aage et al. (2017) "Giga-voxel computational morphogenesis for structural design"

<img src="/_images/mechdesign/featop5.png" alt="Drawing" style="height: 100px;"/> 

Image from Wu et al. (2018) "Infill Optimization for Additive Manufacturing - Approaching Bone-like Porous Structures"

<img src="/_images/mechdesign/featop6.png" alt="Drawing" style="height: 100px;"/> 

Image from Dassault Systems



