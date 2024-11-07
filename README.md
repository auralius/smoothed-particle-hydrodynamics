[![View Smoothed particle hydrodynamics on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/96334-smoothed-particle-hydrodynamics)

# smoothed-particle-hydrodynamics
A MATLAB implementation of the smoothed particle hydrodynamics based on the Philip Mocz's paper, in both CPU and GPU. 

By using a GPU, we can speed up the computation to about 10 times. Here, I am using the NVidia Tesla K20m. It is probably the cheapest Tesla card in second-hand markets. ;-)
A normal GPU also works.  

Please read Philip Mocz's paper here:  
https://pmocz.github.io/manuscripts/pmocz_sph.pdf  

Although the paper stated MATLAB implementation, I can only find the Python implementation in his Github repository:  
https://github.com/pmocz/sph-python


These are some examples in this repository:  

Example 1:  
Free-falling Universitas Telkom's logo  
![alt text](https://github.com/auralius/smoothed-particle-hydrodynamics/blob/main/figures/sph_demo_telu.gif)


Example 2:  
Free-falling particles arranged in a grid formation  
![alt text](https://github.com/auralius/smoothed-particle-hydrodynamics/blob/main/figures/sph_demo1a.gif)

Example 3:  
More free-falling particles also arranged in a grid formation  ;-)  
![alt text](https://github.com/auralius/smoothed-particle-hydrodynamics/blob/main/figures/sph_demo1b.gif)

Example 4:  
Waaayy more free-falling particles also arranged in a grid formation  ;-) ;-)    
![alt text](https://github.com/auralius/smoothed-particle-hydrodynamics/blob/main/figures/sph_demo1c.gif)


Additional notes:

SPHDemo2D_Ex2_NeighbourSearch_CPU.m uses a very simple neighbour search mechanism to reduce the execution time. We first sort all nodes based on their distances to the origin. Hence, the adjacent nodes are neighbours.

Contact:  
manurunga@yadex.com

