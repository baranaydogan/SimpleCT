/* 
------------------------------------------------------------
Dogu Baran Aydogan - baran.aydogan@gmail.com
21.07.2014
------------------------------------------------------------
 
------------Computes and prunes contour tree----------------

- Carr's split and merge algorithm is used
(H. Carr, J. Snoeyink, and U. Axen. Computing contour
trees in all dimensions. In Proceedings of the
eleventh annual symposium on Discrete algorithms, pages
918–926, San Francisco, California, United States,
2000. Society for Industrial and Applied Mathematic.)

- A union-find like data structure is used to store and manipulate
join, split and contour trees

- Pruning is done based on similar method explained in my
doctoral thesis (Dogu Baran Aydogan. Contour Tree Connectivity
and Analysis of Microstructures. Tampere University of
Technology. 2014)

- Pruning simply removes edges whose area is below an area
threshold or if the intensity gap between the vertices of the
edges are below an intensity threshold.

- After join and split trees are merged, an intermediate tree
named the "merge tree (MT)" is computed. This is actually the
non-pruned contour tree. However since it is not pruned it is not
very useful. The prune function operates on MT and
generates the "contour tree (CT)" which is the main output.

- Since prune function does not alter merge tree, it can be called
with different area and intensity thresholds to compute a useful
contour tree.
_____________________________________________________________

EXAMPLE USE:

double* img;
// Fill contents of img

vector<unsigned int> size;
// Fill with dimensions of the image so that
// size[0] = number of rows
// size[1] = number of columns
// size[2] = number of slices (optional)

// The constructor creates the merge tree immediately
ContourTree ct(img,size);

// Prune the merge tree with area threshold 16 and intensity
// threshold 0.1
ct.prune(16, 0.1);

// If desired pruned the merge tree with other parameters
ct.prune(32, 0.4);


_____________________________________________________________
WARNING:
The author does not accept any responsibility or liability for the 
accuracy of the output since it has not been tested exhaustively. 
Please use at your own risk.

_____________________________________________________________
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
_____________________________________________________________

 */

#ifndef SIMPLECT_H_
#define SIMPLECT_H_

#include <utility>
#include <vector>

typedef unsigned int 		uint;
typedef unsigned short 	ushort;
typedef signed char 		schar;

class Element;
class Vertex;
class Leader;
class CT_Vertex;

/*
 * Each element (pixels or voxels) is stored with the
 * following data type
 */

class Element {

public:
	uint 		JT_parent; 	// This element comes from JT_parent in the join tree
	Vertex*		JT_root; 	// This element belongs to JT_root
	uint		ST_parent;	// This element comes from ST_parent in the split tree
	Vertex*		ST_root; 	// This element belongs to ST_root
	uint		MT_parent; 	// This element comes from MT_parent in the merge tree
	uint		CT_parent; 	// This element comes from CT_parent in the contour tree
	CT_Vertex*	CT_root; 	// This element belongs to CT_root

	~Element(){ }
};

/*
 * The vertices of join tree (JT), split tree (ST) and merge tree (MT)
 * are stored with the following data type
 */

class Vertex {

public:
	uint			index; 			// A vertex has the index of the element
	Leader* 		leader; 		// It has a Leader which points to the last element of the tree (this is the union-find like part of the data structure)
	Vertex* 		toVertex; 		// A vertex points to one other vertex - which creates an edge. A vertex cannot point to multiple vertices or no vertex.
	short			childrenCount; 	// A vertex has a number of children. If it is 0 then the edge is a pendant.
	std::vector<Vertex*>	childrenVertices; 	// A vertex has the list of children Vertices.

	~Vertex(){childrenVertices.clear();}
};

/*
 * The vertices of join tree (JT), split tree (ST) and merge tree (MT)
 * has a Leader that points to the last element of the tree.
 * The Leader leads all vertices that have been swept before and
 * connected. Therefore, for the next element to be swept, one does
 * not need to trace all the JT_parent, ST_parent or MT_parent to
 * their last index. (This is the union-find like data structure of this
 * code.)
 */

class Leader {

public:
	uint			index; 			// Index of the element that is lead to
	std::vector<Vertex*>	vertices;		// The list of vertices that this leader leads

	~Leader(){vertices.clear();}
};

/*
 * Contour tree has a its own type of vertex structure
 */


class CT_Vertex {

public:
	uint 				index; 		// Element index of the vertex
	CT_Vertex*			toVertex;	// A contour tree vertex always points to one other contour tree vertex. It cannot point multiple vertices or no vertex.
	uint 				area; 		// Since this vertex goes to one other vertex, it creates an edge. The number of elements on this edge is stored as area.
	std::vector<CT_Vertex*>	fromVertices; 	// The list of vertices that point to this vertex.

	~CT_Vertex(){fromVertices.clear();}
};

/*
 * Contour tree class
 */

class ContourTree {

public:

	ContourTree(double* inputImg, std::vector<uint> inputImgSize); 	// Constructor requires the image pointer and the dimension of the image
	~ContourTree();

	void 		prune(uint areaThresh, double intensityThresh);	// Prunes the merge tree and creates the contour tree
	uint*       export_CT(); 						// Exports the contour tree to the caller (read declaration for detail)
	uint* 		export_CT_img(); 					// Exports an image for the contour tree to the caller. Here the intensity of pixels is equal to the index of the contour tree vertex that the pixel belongs to.

	void 		print_CT(char * ctFileName);				// Prints the contour tree graph on a text file. The graph can be plotted in MATLAB using "plotCT" function.
	void 		print_CT_img(char* ctImgFileName); 			// An image (CT_img.raw) for the contour tree is written in raw format. Here the intensity of pixels is equal to the index of the contour tree vertex that the pixel belongs to.
	void 		print_Short(); 					// Writes the number of merge and contour tree vertices

	double* 	getImg() {return img;}; 				// Returns the input image
	ushort 		getImgDimension() {return dimension;}; 		// Returns the dimension
	uint* 		getImgSize() {return size;};				// Returns the sizes of the dimensions
	uint* 		getSortedImgIndices() {return indices;}; 		// Returns the sorted indices
	std::vector<CT_Vertex*>* getCT() {return &CT;}; 			// Returns the contour tree vector

private:

	void 		perturbAndSortIndices(); 				// The algorithm starts with perturbing and sorting the indices
	void 		compute_JT(); 						// The first sweep computes the join tree
	void 		compute_ST(); 						// The second sweep computes the split tree
	void		augment_JT(); 						// Augments the join tree
	void		augment_ST(); 						// Augments the split tree
	void 		merge(); 						// Augmented trees are merged and the merge tree (MT) is generated

	double*		img;
	ushort 		dimension;
	uint*		size;
	uint*		indices;

	std::vector<Element*> 				elements; 		// Vector of elements
	std::vector<Vertex*> 				JT; 			// The join tree as a vector of vertices
	std::vector<Vertex*> 				ST; 			// The split tree as a vector of vertices
	std::vector< std::pair<Vertex*, bool> >	MT;			// Merge tree of pair data type. The first element is the vertex, the second is logical marking if the vertex is from join or split tree.
	std::vector<CT_Vertex*> CT; 						// Contour tree as a vector of contour tree vertices
};

#endif
