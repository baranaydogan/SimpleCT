SimpleCT
========

SimpleCT is a simple contour tree library with MATLAB interface.

If you use this code in your work please cite the following publication:

Aydogan D.B. & Hyttinen J. “Contour tree connectivity of binary images from 
algebraic graph theory”, IEEE International Conference on Image Processing (ICIP), 
15-22.09.2013, Melbourne, Australia

-------------------------------------------------------------------
DESCRIPTION:

SimpleCT simply computes the contour tree of 2D/3D images and does
heuristic pruning.

Contour tree computation algorithm is based on the paper:

H. Carr, J. Snoeyink, and U. Axen. Computing contour trees in all 
dimensions. In Proceedings of the eleventh annual symposium on 
Discrete algorithms, pages 918–926, San Francisco, California, 
United States, 2000. Society for Industrial and Applied Mathematic.

This implementation uses OpenMP for multicore processing. With this,
join and split trees are computed in parallel. Although, it is not
required, it is strongly recommended to install OpenMP and link the
libraries for faster computation.

-------------------------------------------------------------------
FAST GUIDE:

If you are familiar with MATLAB or intend to use SimpleCT in MATLAB,
please go to the "matlab" folder and read the README file. This 
explains and directs you on how to install and use SimpleCT.

If you are planning to use SimpleCT in C/C++, please go to "example"
folder and read the README file.

____________________________________________________________________
LICENSE INFORMATION:

SimpleCT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by 
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimpleCT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
____________________________________________________________________

Dogu Baran Aydogan - baran.aydogan@gmail.com
22.07.2014
____________________________________________________________________
