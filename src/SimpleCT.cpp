/*
 * _____________________________________________________________________
 * SimpleCT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SimpleCT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * _____________________________________________________________________
 */

#include <stddef.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <climits>
#include <cmath>

#include "SimpleCT.h"

bool sortByIntensity ( const std::pair<double,uint> &l, const std::pair<double,uint> &r) { return l.first < r.first; }

void ContourTree::perturbAndSortIndices() {
    
    uint N = 1;
    for (uint i=0; i<dimension; ++i)
        N = N*size[i];
    
    std::pair<double,uint> *out;
    out 		= new std::pair<double,uint> [N];
    
    for (uint i=0; i<N; ++i)
        out[i] = std::make_pair( img[i]+0.000001*i/(N-1), i );
    
    for ( uint i = 0; i < N; ++i )
        img[i] = out[i].first;
    
    std::sort ( out, out+N, sortByIntensity);
    
    for ( uint i = 0; i < N; ++i )
        indices[i] = out[i].second;
    
    delete[] out;
    
}

void getNeighborIndices(std::vector<uint> *neighborIndices, uint ind, uint* size, ushort dimension) {
    
    if (dimension > 2) {
        
        int m = size[0];
        int n = size[1];
        int p = size[2];
        
        int x,	y,	z, 	mm, nm, pm;
        int xm,	ym,	zm, xp, yp, zp;
        
        mm = m - 1;
        nm = n - 1;
        pm = p - 1;
        
        z 	=	ind / (m*n);
        ind	%= 	(m*n);
        y 	=	ind / m;
        x 	=	ind % m;
        
        xm 	= x - 1;
        xp 	= x + 1;
        ym 	= y - 1;
        yp 	= y + 1;
        zm 	= z - 1;
        zp 	= z + 1;
        
        
        if ( x < mm )							neighborIndices->push_back( xp+y*m+z*m*n );
        if ( x > 0 )        					neighborIndices->push_back( xm+y*m+z*m*n );
        if ( y < nm )		   					neighborIndices->push_back( x+yp*m+z*m*n );
        if ( y > 0 )           					neighborIndices->push_back( x+ym*m+z*m*n );
        if ( z < pm )   						neighborIndices->push_back( x+y*m+zp*m*n );
        if ( z > 0 )           					neighborIndices->push_back( x+y*m+zm*m*n );
        
        if ( x < mm && y < nm )					neighborIndices->push_back( xp+yp*m+z*m*n );
        if ( x > 0	&& y < nm )					neighborIndices->push_back( xm+yp*m+z*m*n );
        if ( x < mm && y > 0  )					neighborIndices->push_back( xp+ym*m+z*m*n );
        if ( x > 0	&& y > 0  )					neighborIndices->push_back( xm+ym*m+z*m*n );
        
        if ( x < mm && z < pm )					neighborIndices->push_back( xp+y*m+zp*m*n );
        if ( x > 0 	&& z < pm )					neighborIndices->push_back( xm+y*m+zp*m*n );
        if ( x < mm && z > 0  )					neighborIndices->push_back( xp+y*m+zm*m*n );
        if ( x > 0 	&& z > 0  )					neighborIndices->push_back( xm+y*m+zm*m*n );
        
        if ( y < nm	&& z < pm )					neighborIndices->push_back( x+yp*m+zp*m*n );
        if ( y > 0	&& z < pm )					neighborIndices->push_back( x+ym*m+zp*m*n );
        if ( y < nm	&& z > 0  )					neighborIndices->push_back( x+yp*m+zm*m*n );
        if ( y > 0	&& z > 0  ) 				neighborIndices->push_back( x+ym*m+zm*m*n );
        
        if ( x < mm	&& y < nm && z < pm )		neighborIndices->push_back( xp+yp*m+zp*m*n );
        if ( x > 0	&& y < nm && z < pm )		neighborIndices->push_back( xm+yp*m+zp*m*n );
        if ( x < mm	&& y > 0  && z < pm )		neighborIndices->push_back( xp+ym*m+zp*m*n );
        if ( x > 0	&& y > 0  && z < pm )		neighborIndices->push_back( xm+ym*m+zp*m*n );
        
        if ( x < mm && y < nm && z > 0 )		neighborIndices->push_back( xp+yp*m+zm*m*n );
        if ( x > 0 	&& y < nm && z > 0 )		neighborIndices->push_back( xm+yp*m+zm*m*n );
        if ( x < mm && y > 0  && z > 0 )		neighborIndices->push_back( xp+ym*m+zm*m*n );
        if ( x > 0	&& y > 0  && z > 0 )		neighborIndices->push_back( xm+ym*m+zm*m*n );
        
    } else {
        
        int m = size[0];
        int n = size[1];
        
        int x,	y,	mm, nm;
        int xm,	ym, xp, yp;
        
        mm = m - 1;
        nm = n - 1;
        
        y 	=	ind / m;
        x 	=	ind % m;
        
        xm 	= x - 1;
        xp 	= x + 1;
        ym 	= y - 1;
        yp 	= y + 1;
        
        if ( x < mm )							neighborIndices->push_back( xp+y*m );
        if ( x > 0 	)        					neighborIndices->push_back( xm+y*m );
        if ( y < nm )   						neighborIndices->push_back( x+yp*m );
        if ( y > 0 	)           				neighborIndices->push_back( x+ym*m );
        
        if ( x < mm	&& y < nm 	)				neighborIndices->push_back( xp+yp*m );
        if ( x > 0	&& y > 0 	)				neighborIndices->push_back( xm+ym*m );
        if ( x < mm	&& y > 0 	)				neighborIndices->push_back( xp+ym*m );
        if ( x > 0	&& y < nm 	)				neighborIndices->push_back( xm+yp*m );
        
    }
    
    
}

void ContourTree::compute_JT(){
    
    uint N = 1;
    for (uint i=0; i<dimension; ++i)
        N = N*size[i];
    
    std::vector<uint> 			neighborIndices;
    std::vector<uint>::iterator neighborIndices_it;
    
    std::set<uint> 				mergingRegions;
    std::set<uint>::iterator 	mergingRegions_it;
    
    Leader* tmpLeader;
    Vertex* tmpVertex;
    Leader* deleteLeader;
    
    std::vector<Vertex*>::iterator vertex_it;
    
    for (uint i = N-1; i<UINT_MAX; --i) {
        
        getNeighborIndices(&neighborIndices, indices[i], size, dimension);
        for (neighborIndices_it = neighborIndices.begin(); neighborIndices_it != neighborIndices.end(); ++neighborIndices_it)
            if (elements[*neighborIndices_it]->JT_root)
                mergingRegions.insert(elements[*neighborIndices_it]->JT_root->leader->index);
        
        if (mergingRegions.size() == 1) {
            elements[*mergingRegions.begin()]->JT_parent = indices[i];
            elements[indices[i]]->JT_root = elements[*mergingRegions.begin()]->JT_root;
            elements[*mergingRegions.begin()]->JT_root->leader->index = indices[i];
        } else if ( mergingRegions.empty() ) {
            tmpLeader						= new Leader;
            tmpVertex						= new Vertex;
            tmpLeader->index 				= indices[i];
            tmpVertex->leader 				= tmpLeader;
            tmpVertex->index 				= indices[i];
            tmpVertex->childrenCount 		= 0;
            tmpLeader->vertices.push_back(tmpVertex);
            elements[indices[i]]->JT_root 	= tmpVertex;
            JT.push_back(tmpVertex);
        }else {
            tmpLeader						= new Leader;
            tmpVertex						= new Vertex;
            tmpLeader->index 				= indices[i];
            tmpVertex->leader 				= tmpLeader;
            tmpVertex->index 				= indices[i];
            tmpVertex->childrenCount 		= 0;
            tmpLeader->vertices.push_back(tmpVertex);
            elements[indices[i]]->JT_root 	= tmpVertex;
            JT.push_back(tmpVertex);
            
            tmpVertex->toVertex 			= NULL;
            for (mergingRegions_it = mergingRegions.begin(); mergingRegions_it != mergingRegions.end(); ++mergingRegions_it) {
                elements[*mergingRegions_it]->JT_parent = indices[i];
                tmpVertex->childrenVertices.push_back(elements[*mergingRegions_it]->JT_root);
                (tmpVertex->childrenCount)++;
                deleteLeader				= elements[*mergingRegions_it]->JT_root->leader;
                tmpLeader->vertices.insert(tmpLeader->vertices.end(),deleteLeader->vertices.begin(),deleteLeader->vertices.end() );
                for (vertex_it = deleteLeader->vertices.begin(); vertex_it != deleteLeader->vertices.end(); ++vertex_it)
                    (*vertex_it)->leader 	= tmpLeader;
                delete deleteLeader;
                elements[*mergingRegions_it]->JT_root->toVertex = tmpVertex;
            }
        }
        
        mergingRegions.clear();
        neighborIndices.clear();
        
    }
    
}

void ContourTree::compute_ST(){
    
    uint N = 1;
    for (uint i=0; i<dimension; ++i)
        N = N*size[i];
    
    std::vector<uint> 			neighborIndices;
    std::vector<uint>::iterator neighborIndices_it;
    
    std::set<uint> 				mergingRegions;
    std::set<uint>::iterator 	mergingRegions_it;
    
    Leader* tmpLeader;
    Vertex* tmpVertex;
    Leader* deleteLeader;
    
    std::vector<Vertex*>::iterator vertex_it;
    
    for (uint i = 0; i < N; ++i) {
        
        getNeighborIndices(&neighborIndices, indices[i], size, dimension);
        
        for (neighborIndices_it = neighborIndices.begin(); neighborIndices_it != neighborIndices.end(); ++neighborIndices_it)
            if (elements[*neighborIndices_it]->ST_root)
                mergingRegions.insert(elements[*neighborIndices_it]->ST_root->leader->index);
        
        if (mergingRegions.size() == 1) {
            elements[*mergingRegions.begin()]->ST_parent = indices[i];
            elements[indices[i]]->ST_root = elements[*mergingRegions.begin()]->ST_root;
            elements[*mergingRegions.begin()]->ST_root->leader->index = indices[i];
        } else if ( mergingRegions.empty() ) {
            tmpLeader						= new Leader;
            tmpVertex						= new Vertex;
            tmpLeader->index 				= indices[i];
            tmpVertex->leader 				= tmpLeader;
            tmpVertex->index 				= indices[i];
            tmpVertex->childrenCount 		= 0;
            tmpLeader->vertices.push_back(tmpVertex);
            elements[indices[i]]->ST_root 	= tmpVertex;
            ST.push_back(tmpVertex);
        }else {
            tmpLeader						= new Leader;
            tmpVertex						= new Vertex;
            tmpLeader->index 				= indices[i];
            tmpVertex->leader 				= tmpLeader;
            tmpVertex->index 				= indices[i];
            tmpVertex->childrenCount 		= 0;
            tmpLeader->vertices.push_back(tmpVertex);
            elements[indices[i]]->ST_root 	= tmpVertex;
            tmpVertex->toVertex 			= NULL;
            ST.push_back(tmpVertex);
            for (mergingRegions_it = mergingRegions.begin(); mergingRegions_it != mergingRegions.end(); ++mergingRegions_it) {
                elements[*mergingRegions_it]->ST_parent = indices[i];
                tmpVertex->childrenVertices.push_back(elements[*mergingRegions_it]->ST_root);
                (tmpVertex->childrenCount)++;
                deleteLeader				= elements[*mergingRegions_it]->ST_root->leader;
                tmpLeader->vertices.insert(tmpLeader->vertices.end(),deleteLeader->vertices.begin(),deleteLeader->vertices.end() );
                for (vertex_it = deleteLeader->vertices.begin(); vertex_it != deleteLeader->vertices.end(); ++vertex_it)
                    (*vertex_it)->leader 	= tmpLeader;
                delete deleteLeader;
                elements[*mergingRegions_it]->ST_root->toVertex = elements[indices[i]]->ST_root;
            }
        }
        
        mergingRegions.clear();
        neighborIndices.clear();
        
    }
    
}

ContourTree::ContourTree(double* inputImg, std::vector<uint>  inputImgSize) {
    
    dimension 	= inputImgSize.size();
    size 		= new uint[dimension];
    
    uint numberOfElements = 1;
    
    for (ushort i = 0; i < dimension; ++i) {
        size[i] 			= inputImgSize[i];
        numberOfElements 	= numberOfElements*inputImgSize[i];
    }
    
    indices 	= new uint[numberOfElements];
    
    img 		= new double[numberOfElements];
    
    Element* tmp;
    
    for(uint i = 0; i < numberOfElements; ++i) {
        tmp 			= new Element();
        tmp->JT_parent 	= UINT_MAX;
        tmp->JT_root	= NULL;
        tmp->ST_parent 	= UINT_MAX;
        tmp->ST_root	= NULL;
        
        tmp->MT_parent 	= UINT_MAX;
        elements.push_back(tmp);
        
        img[i] 			= double(inputImg[i]);
    }
    
    perturbAndSortIndices();
    
#pragma omp parallel sections
{
#pragma omp section
{
    compute_JT();
}
#pragma omp section
{
    compute_ST();
}
}

this->augment_JT();
this->augment_ST();

this->merge();

}

ContourTree::~ContourTree() {
    
    if (JT.size() != 0)
        delete JT.back()->leader;
    
    if (ST.size() != 0)
        delete ST.back()->leader;
    
    std::vector<Vertex*>::iterator itV;
    for(itV = JT.begin(); itV != JT.end(); ++itV)
        delete *itV;
    JT.clear();
    
    for(itV = ST.begin(); itV != ST.end(); ++itV)
        delete *itV;
    ST.clear();
    
    std::vector<Element*>::iterator itE;
    for(itE = elements.begin(); itE != elements.end(); ++itE)
        delete *itE;
    elements.clear();
    
    std::vector<CT_Vertex*>::iterator itCT;
    for(itCT = CT.begin(); itCT != CT.end(); ++itCT)
        delete *itCT;
    CT.clear();
    
    delete[] size;
    delete[] indices;
    delete[] img;
    
}

void ContourTree::augment_JT() {
    
    Vertex* fromVertex;
    Vertex* tmp;
    Element* tmpEl;
    
    std::vector<Vertex*>::iterator it;
    
    for (it = ST.begin(); it < ST.end(); ++it)
        if (elements[(*it)->index]->JT_root->index != (*it)->index) {
        
        tmpEl 				= elements[(*it)->index];
        
        fromVertex 			= tmpEl->JT_root;
        tmp 				= new Vertex;
        
        tmp->index 			= (*it)->index;
        tmp->toVertex 		= fromVertex->toVertex;
        tmp->leader 		= fromVertex->leader;
        tmp->childrenCount 	= 1;
        
        tmp->childrenVertices.push_back(fromVertex);
        
        tmpEl->JT_root 		= tmp;
        
        if (tmpEl->JT_parent == UINT_MAX)
            tmpEl->JT_parent = (*it)->index;
        
        if (tmp->toVertex) {
            
            while (tmpEl->JT_parent != tmp->toVertex->index) {
                tmpEl 			= elements[tmpEl->JT_parent];
                tmpEl->JT_root 	= tmp;
            }
            
            for(uint i=0; i<fromVertex->toVertex->childrenVertices.size();++i)
                if(fromVertex->toVertex->childrenVertices.at(i) == fromVertex) {
                fromVertex->toVertex->childrenVertices.erase(fromVertex->toVertex->childrenVertices.begin()+i);
                break;
                }
            
            fromVertex->toVertex->childrenVertices.push_back(tmp);
            
        } else {
            
            while (tmpEl->JT_parent != tmp->leader->index) {
                tmpEl 			= elements[tmpEl->JT_parent];
                tmpEl->JT_root 	= tmp;
            }
            tmpEl 			= elements[tmpEl->JT_parent];
            tmpEl->JT_root 	= tmp;
            
        }
        
        fromVertex->toVertex = tmp;
        JT.push_back(tmp);
        tmp->leader->vertices.push_back(tmp);
        
        }
    
}

void ContourTree::augment_ST() {
    
    Vertex* fromVertex;
    Vertex* tmp;
    Element* tmpEl;
    
    std::vector<Vertex*>::iterator it;
    
    for (it = JT.begin(); it < JT.end(); ++it)
        if (elements[(*it)->index]->ST_root->index != (*it)->index) {
        
        tmpEl 				= elements[(*it)->index];
        
        fromVertex 			= tmpEl->ST_root;;
        tmp 				= new Vertex;
        
        tmp->index 			= (*it)->index;
        tmp->toVertex 		= fromVertex->toVertex;
        tmp->leader 		= fromVertex->leader;
        tmp->childrenCount 	= 1;
        
        tmp->childrenVertices.push_back(fromVertex);
        
        tmpEl->ST_root 		= tmp;
        
        if (tmpEl->ST_parent == UINT_MAX)
            tmpEl->ST_parent = (*it)->index;
        
        if (tmp->toVertex) {
            
            while (tmpEl->ST_parent != tmp->toVertex->index) {
                tmpEl 			= elements[tmpEl->ST_parent];
                tmpEl->ST_root 	= tmp;
            }
            
            for(uint i=0; i<fromVertex->toVertex->childrenVertices.size();++i)
                if(fromVertex->toVertex->childrenVertices.at(i) == fromVertex) {
                fromVertex->toVertex->childrenVertices.erase(fromVertex->toVertex->childrenVertices.begin()+i);
                break;
                }
            
            fromVertex->toVertex->childrenVertices.push_back(tmp);
            
            
        } else {
            
            while (tmpEl->ST_parent != tmp->leader->index) {
                tmpEl 			= elements[tmpEl->ST_parent];
                tmpEl->JT_root 	= tmp;
            }
            tmpEl 			= elements[tmpEl->ST_parent];
            tmpEl->ST_root 	= tmp;
            
        }
        
        fromVertex->toVertex = tmp;
        ST.push_back(tmp);
        tmp->leader->vertices.push_back(tmp);
        
        }
    
}

void ContourTree::merge() {
    
    bool stopFlag 		= false;
    
    while(!stopFlag) {
        
        stopFlag 		= true;
        
#pragma omp parallel for shared(stopFlag)
for (uint i=0; i<JT.size(); i++)
    if ( ( JT.at(i)->childrenCount == 0 ) && ( elements[JT.at(i)->index]->ST_root->childrenCount == 1 ) && ( JT.at(i)->toVertex != NULL ) ) {
    
    while(elements[JT.at(i)->toVertex->index]->MT_parent != UINT_MAX)
        JT.at(i)->toVertex = JT.at(i)->toVertex->toVertex;
    
    
    Element* tmp_El = elements[JT.at(i)->index];
    Element* tmp_it_El;
    
    while (tmp_El != elements[JT.at(i)->toVertex->index])
        if (elements[tmp_El->JT_parent]->MT_parent == UINT_MAX) {
        tmp_El->MT_parent	= tmp_El->JT_parent;
        tmp_El 					= elements[tmp_El->JT_parent];
        }
        else {
        tmp_it_El = tmp_El;
        while ( (elements[tmp_it_El->JT_parent]->MT_parent != UINT_MAX) && (tmp_it_El->JT_parent != JT.at(i)->toVertex->index) )
            tmp_it_El = elements[tmp_it_El->JT_parent];
        tmp_El->MT_parent	= tmp_it_El->JT_parent;
        tmp_El 					= elements[tmp_it_El->JT_parent];
        }
    
#pragma omp critical
MT.push_back( std::make_pair(JT.at(i),true) );
JT.at(i)->childrenCount 					= -1;
(JT.at(i)->toVertex->childrenCount)--;
stopFlag 									= false;
    }

#pragma omp parallel for shared(stopFlag)
for (uint i=0; i<ST.size(); i++)
    if ( ( ST.at(i)->childrenCount == 0 ) && ( elements[ST.at(i)->index]->JT_root->childrenCount == 1 ) && ( ST.at(i)->toVertex != NULL ) ) {
    
    while(elements[ST.at(i)->toVertex->index]->MT_parent != UINT_MAX)
        ST.at(i)->toVertex = ST.at(i)->toVertex->toVertex;
    
    Element* tmp_El = elements[ST.at(i)->index];
    Element* tmp_it_El;
    
    while (tmp_El != elements[ST.at(i)->toVertex->index])
        if (elements[tmp_El->ST_parent]->MT_parent == UINT_MAX) {
        tmp_El->MT_parent	= tmp_El->ST_parent;
        tmp_El 					= elements[tmp_El->ST_parent];
        }
        else {
        
        tmp_it_El = tmp_El;
        
        while ( (elements[tmp_it_El->ST_parent]->MT_parent != UINT_MAX) && (tmp_it_El->ST_parent != ST.at(i)->toVertex->index) )
            tmp_it_El = elements[tmp_it_El->ST_parent];
        
        tmp_El->MT_parent	= tmp_it_El->ST_parent;
        tmp_El 					= elements[tmp_it_El->ST_parent];
        }
    
#pragma omp critical
MT.push_back( std::make_pair(ST.at(i),false) );
ST.at(i)->childrenCount 					= -1;
(ST.at(i)->toVertex->childrenCount)--;
stopFlag 									= false;
    }

    }
    
    
}

void ContourTree::prune(uint areaThresh, double intensityThresh) {
    
    for(uint i = 0; i<MT.size(); i++){
        
        CT_Vertex* tmp_CT_Vertex 	= new CT_Vertex;
        
        CT.push_back(tmp_CT_Vertex);
        
        tmp_CT_Vertex->index 		= MT.at(i).first->index;
        
        Element* tmp_Element 		= elements[MT.at(i).first->index];
        uint area 					= 1;
        
        
        tmp_Element->CT_parent 	= tmp_Element->MT_parent;
        tmp_Element->CT_root 	= tmp_CT_Vertex;
        
        while(tmp_Element->MT_parent 	!= MT.at(i).first->toVertex->index){
            tmp_Element 				= elements[tmp_Element->MT_parent];
            tmp_Element->CT_parent 		= tmp_Element->MT_parent;
            tmp_Element->CT_root 		= tmp_CT_Vertex;
            area++;
        }
        
        tmp_CT_Vertex->area	= area;
    }
    
    uint numberOfVertices = CT.size();
    
    for(uint i = 0; i<numberOfVertices; i++){
        
        uint index 				= MT.at(i).first->toVertex->index;
        bool vertexFound 		= false;
        
        for(uint j = 0; j<CT.size(); j++)
            if (CT.at(j)->index == index) {
            CT.at(i)->toVertex 	= CT.at(j);
            vertexFound 		= true;
            break;
            }
        
        if(vertexFound == false) {
            CT_Vertex* tmp_CT_Vertex 	= new CT_Vertex;
            CT.push_back(tmp_CT_Vertex);
            tmp_CT_Vertex->index 		= index;
            elements[index]->CT_parent  = index;
            elements[index]->CT_root    = tmp_CT_Vertex;
            tmp_CT_Vertex->toVertex 	= tmp_CT_Vertex;
            tmp_CT_Vertex->area 		= 1;
            CT.at(i)->toVertex 			= tmp_CT_Vertex;
        }
    }
    
    for(uint i = 0; i<CT.size(); i++){
        uint index = CT.at(i)->index;
        for(uint j = 0; j<CT.size(); j++)
            if (CT.at(j)->toVertex->index == index)
                CT.at(i)->fromVertices.push_back(CT.at(j));
    }
    
    
    uint N = 1;
    for (uint i=0; i<dimension; ++i)
        N = N*size[i];
    
    std::vector<CT_Vertex*> tmp_CT_vec;
    std::vector<CT_Vertex*>::iterator CT_it;
    std::vector<CT_Vertex*>::iterator tmp_CT_it;
    std::vector<CT_Vertex*>::iterator CT_children_it;
    
    CT_Vertex* tmp_CT_fromVertex;
    CT_Vertex* tmp_CT_toVertex;
    
    Element* tmp_El;
    
    uint prev_ind;
    uint cur_ind;
    uint next_ind;
    
    double fromIntensity;
    double toIntensity;
    
    double intensityDiff;
    bool diffSign;
    
    double comingRefIntensity;
    
    CT_Vertex* appendToComing;
    CT_Vertex* appendToGoing;
    
    
    bool contIteration = true;
    
    while (contIteration) {
        
        uint numberOfVerticesRemoved 	= 1;
        bool firstStageRemovalDone 		= false;
        
        while (numberOfVerticesRemoved > 0) {
            
            numberOfVerticesRemoved = 0;
            
            for(CT_it = CT.begin(); CT_it != CT.end(); CT_it++) {
                
                bool thisVertexRemoved = false;
                
                if ( (*CT_it)->fromVertices.size() == 0) {
                    
                    if ( (*CT_it)->toVertex->fromVertices.size() == 1 ) {
                        
                        // Vertex collapse
                        
                        if ((!thisVertexRemoved) && ( (*CT_it)->toVertex->toVertex != (*CT_it)->toVertex ) ) {
                            
                            tmp_CT_toVertex 	= (*CT_it)->toVertex->toVertex;
                            tmp_CT_fromVertex 	= (*CT_it);
                            
                            tmp_El 				= elements[(*CT_it)->toVertex->index];
                            
                            do {
                                tmp_El->CT_root = tmp_CT_fromVertex;
                                tmp_El = elements[tmp_El->CT_parent];
                            } while( tmp_El->CT_root != tmp_CT_toVertex);
                            
                            std::vector<CT_Vertex*>::iterator tmp_CT_vertex = tmp_CT_toVertex->fromVertices.begin();
                            while( (*tmp_CT_vertex) != tmp_CT_fromVertex->toVertex )
                                tmp_CT_vertex++;
                            tmp_CT_toVertex->fromVertices.erase(tmp_CT_vertex);
                            
                            tmp_CT_toVertex->fromVertices.push_back(tmp_CT_fromVertex);
                            tmp_CT_fromVertex->area = tmp_CT_fromVertex->area + (*CT_it)->toVertex->area;
                            
                            tmp_CT_vertex = CT.begin();
                            while( (*tmp_CT_vertex) != tmp_CT_fromVertex->toVertex )
                                tmp_CT_vertex++;
                            CT.erase(tmp_CT_vertex);
                            
                            delete (*CT_it)->toVertex;
                            
                            tmp_CT_fromVertex->toVertex = tmp_CT_toVertex;
                            
                            numberOfVerticesRemoved++;
                            firstStageRemovalDone = true;
                            thisVertexRemoved = true;
                        }
                        
                    } else {
                        
                        // Checking pendant vertices
                        
                        tmp_CT_fromVertex 	= *CT_it;
                        tmp_CT_toVertex 	= (*CT_it)->toVertex;
                        
                        fromIntensity 		= img[tmp_CT_fromVertex->index];
                        toIntensity 		= img[tmp_CT_toVertex->index];
                        
                        intensityDiff 		= fromIntensity - toIntensity;
                        diffSign 			= intensityDiff > 0;
                        
                        if ( (!thisVertexRemoved) && ( (tmp_CT_fromVertex->index != indices[N-1]) && (tmp_CT_fromVertex->index != indices[0]) ) && (tmp_CT_fromVertex->area <= areaThresh || std::abs(intensityDiff) <= intensityThresh) ) {
                            
                            // Check if can be append to going
                            
                            appendToGoing = tmp_CT_toVertex;
                            
                            if ( (diffSign && ( img[tmp_CT_toVertex->toVertex->index] > fromIntensity) ) || (!diffSign && ( img[tmp_CT_toVertex->toVertex->index] < fromIntensity) ) )
                                appendToGoing = tmp_CT_toVertex->toVertex;
                            
                            if (appendToGoing != tmp_CT_toVertex) {
                                
                                tmp_El = elements[tmp_CT_fromVertex->index];
                                
                                while( tmp_El->CT_parent != tmp_CT_toVertex->index) {
                                    tmp_El->CT_root = tmp_CT_toVertex;
                                    tmp_El = elements[tmp_El->CT_parent];
                                }
                                
                                tmp_El->CT_root = tmp_CT_toVertex;
                                tmp_El->CT_parent = elements[tmp_CT_toVertex->index]->CT_parent;
                                
                                elements[tmp_CT_toVertex->index]->CT_parent = tmp_CT_fromVertex->index;
                                
                                CT_children_it = tmp_CT_toVertex->fromVertices.begin();
                                while( (*CT_children_it) != tmp_CT_fromVertex )
                                    CT_children_it++;
                                tmp_CT_toVertex->fromVertices.erase(CT_children_it);
                                
                                tmp_CT_toVertex->area = tmp_CT_toVertex->area + tmp_CT_fromVertex->area;
                                
                                delete *CT_it;
                                CT.erase(CT_it);
                                CT_it--;
                                
                                numberOfVerticesRemoved++;
                                firstStageRemovalDone = true;
                                thisVertexRemoved = true;
                                
                            }
                            
                            
                            // Check if can be append to coming
                            
                            comingRefIntensity 	= fromIntensity;
                            appendToComing 		= tmp_CT_toVertex;
                            
                            for(tmp_CT_it = tmp_CT_toVertex->fromVertices.begin(); tmp_CT_it != tmp_CT_toVertex->fromVertices.end(); tmp_CT_it++)
                                if ( (diffSign && ( img[(*tmp_CT_it)->index] > comingRefIntensity) ) || (!diffSign && ( img[(*tmp_CT_it)->index] < comingRefIntensity) ) ) {
                                comingRefIntensity = img[(*tmp_CT_it)->index];
                                appendToComing = *tmp_CT_it;
                                }
                            
                            if (appendToComing != tmp_CT_toVertex) {
                                
                                tmp_El = elements[tmp_CT_fromVertex->index];
                                
                                while( tmp_El->CT_root != tmp_CT_toVertex) {
                                    tmp_El->CT_root = appendToComing;
                                    tmp_El = elements[tmp_El->CT_parent];
                                }
                                
                                tmp_El = elements[appendToComing->index];
                                
                                while( tmp_El->CT_parent != tmp_CT_toVertex->index)
                                    tmp_El = elements[tmp_El->CT_parent];
                                tmp_El->CT_parent = tmp_CT_fromVertex->index;
                                
                                CT_children_it = tmp_CT_toVertex->fromVertices.begin();
                                while( (*CT_children_it) != tmp_CT_fromVertex )
                                    CT_children_it++;
                                tmp_CT_toVertex->fromVertices.erase(CT_children_it);
                                
                                appendToComing->area = appendToComing->area + tmp_CT_fromVertex->area;
                                
                                delete *CT_it;
                                CT.erase(CT_it);
                                CT_it--;
                                
                                numberOfVerticesRemoved++;
                                firstStageRemovalDone = true;
                                thisVertexRemoved = true;
                                
                            }
                            
                            
                            
                            
                            
                        }
                        
                        
                        if ( (!thisVertexRemoved) && ( (tmp_CT_fromVertex->index != indices[N-1]) && (tmp_CT_fromVertex->index != indices[0]) )&& ( (tmp_CT_fromVertex->area <= areaThresh) && (std::abs(intensityDiff) <= 4*intensityThresh) ) ) {
                            
                            if (tmp_CT_fromVertex->toVertex != tmp_CT_fromVertex) {
                                
                                // Force pendant remove to going edge
                                
                                tmp_El = elements[tmp_CT_fromVertex->index];
                                
                                while( tmp_El->CT_parent != tmp_CT_toVertex->index) {
                                    tmp_El->CT_root = tmp_CT_toVertex;
                                    tmp_El = elements[tmp_El->CT_parent];
                                }
                                
                                tmp_El->CT_root = tmp_CT_toVertex;
                                tmp_El->CT_parent = elements[tmp_CT_toVertex->index]->CT_parent;
                                
                                elements[tmp_CT_toVertex->index]->CT_parent = tmp_CT_fromVertex->index;
                                
                                CT_children_it = tmp_CT_toVertex->fromVertices.begin();
                                while( (*CT_children_it) != tmp_CT_fromVertex )
                                    CT_children_it++;
                                tmp_CT_toVertex->fromVertices.erase(CT_children_it);
                                
                                tmp_CT_toVertex->area = tmp_CT_toVertex->area + tmp_CT_fromVertex->area;
                                
                                delete *CT_it;
                                CT.erase(CT_it);
                                CT_it--;
                                
                                numberOfVerticesRemoved++;
                                firstStageRemovalDone = true;
                                thisVertexRemoved = true;
                                
                            }
                            
                        }
                        
                        
                        
                        
                    }
                    
                }
            }
        }
        
        uint numberOfSaddlesRemoved = 1;
        bool secondStageRemovalDone = false;
        
        while (numberOfSaddlesRemoved > 0) {
            
            numberOfSaddlesRemoved = 0;
            
            for(CT_it = CT.begin(); CT_it != CT.end(); CT_it++) {
                
                bool thisVertexRemoved 		= false;
                
                if ( (!thisVertexRemoved) && ((*CT_it)->fromVertices.size() > 1) && ( (*CT_it)->toVertex != (*CT_it) )  ) {
                    
                    // Saddle collapse
                    
                    tmp_CT_fromVertex 	= *CT_it;
                    tmp_CT_toVertex 	= (*CT_it)->toVertex;
                    
                    fromIntensity 		= img[tmp_CT_fromVertex->index];
                    toIntensity 		= img[tmp_CT_toVertex->index];
                    
                    intensityDiff 		= fromIntensity - toIntensity;
                    
                    if (tmp_CT_fromVertex->area <= areaThresh || std::abs(intensityDiff) <= intensityThresh){
                        
                        tmp_CT_it = tmp_CT_toVertex->fromVertices.begin();
                        while( (*tmp_CT_it) != tmp_CT_fromVertex )
                            tmp_CT_it++;
                        tmp_CT_toVertex->fromVertices.erase(tmp_CT_it);
                        
                        appendToComing 	= tmp_CT_fromVertex->fromVertices.front();
                        
                        for(tmp_CT_it = tmp_CT_fromVertex->fromVertices.begin(); tmp_CT_it != tmp_CT_fromVertex->fromVertices.end(); tmp_CT_it++) {
                            (*tmp_CT_it)->toVertex = tmp_CT_toVertex;
                            
                            tmp_El = elements[(*tmp_CT_it)->index];
                            
                            while(tmp_El->CT_parent != tmp_CT_fromVertex->index)
                                tmp_El = elements[tmp_El->CT_parent];
                            tmp_El->CT_parent = tmp_CT_toVertex->index;
                            
                            tmp_CT_toVertex->fromVertices.push_back(*tmp_CT_it);
                            
                            if ((*tmp_CT_it)->area > appendToComing->area)
                                appendToComing = *tmp_CT_it;
                        }
                        
                        appendToComing->area = appendToComing->area + tmp_CT_fromVertex->area;
                        
                        tmp_El = elements[tmp_CT_fromVertex->index];
                        while(tmp_El != elements[tmp_CT_toVertex->index]) {
                            tmp_El->CT_root 	= appendToComing;
                            tmp_El 				= elements[tmp_El->CT_parent];
                        }
                        
                        tmp_El = elements[appendToComing->index];
                        while(tmp_El->CT_parent != tmp_CT_toVertex->index)
                            tmp_El = elements[tmp_El->CT_parent];
                        tmp_El->CT_parent = tmp_CT_fromVertex->index;
                        
                        delete *CT_it;
                        CT.erase(CT_it);
                        CT_it--;
                        
                        numberOfSaddlesRemoved++;
                        secondStageRemovalDone = true;
                        thisVertexRemoved = true;
                        
                        
                    }
                    
                }
                
                if ( (!thisVertexRemoved) && ( (*CT_it)->fromVertices.size() == 3) && ( (*CT_it)->toVertex->index == (*CT_it)->index ) ) {
                    
                    // Vertex collapse
                    
                    if ((*CT_it)->fromVertices.at(0) == *CT_it) {
                        tmp_CT_toVertex = (*CT_it)->fromVertices.at(1);
                        tmp_CT_fromVertex = (*CT_it)->fromVertices.at(2);
                    } else if ((*CT_it)->fromVertices.at(1) == *CT_it) {
                        tmp_CT_toVertex = (*CT_it)->fromVertices.at(0);
                        tmp_CT_fromVertex = (*CT_it)->fromVertices.at(2);
                    } else {
                        tmp_CT_toVertex = (*CT_it)->fromVertices.at(0);
                        tmp_CT_fromVertex = (*CT_it)->fromVertices.at(1);
                    }
                    
                    if ( ( img[tmp_CT_toVertex->index] > img[(*CT_it)->index] && img[tmp_CT_fromVertex->index] < img[(*CT_it)->index] ) || ( img[tmp_CT_fromVertex->index] > img[(*CT_it)->index] && img[tmp_CT_toVertex->index] < img[(*CT_it)->index] ) ) {
                        
                        prev_ind = tmp_CT_toVertex->index;
                        cur_ind = elements[prev_ind]->CT_parent;
                        
                        do {
                            
                            next_ind = elements[cur_ind]->CT_parent;
                            
                            elements[cur_ind]->CT_parent = prev_ind;
                            elements[cur_ind]->CT_root = tmp_CT_fromVertex;
                            
                            prev_ind = cur_ind;
                            cur_ind = next_ind;
                            
                        } while( prev_ind != (*CT_it)->index );
                        
                        tmp_CT_toVertex->fromVertices.push_back(tmp_CT_fromVertex);
                        tmp_CT_fromVertex->area = tmp_CT_fromVertex->area + tmp_CT_toVertex->area+1;
                        
                        tmp_CT_fromVertex->toVertex = tmp_CT_toVertex;
                        tmp_CT_toVertex->toVertex = tmp_CT_toVertex;
                        tmp_CT_toVertex->fromVertices.push_back(tmp_CT_toVertex);
                        tmp_CT_toVertex->area = 1;
                        
                        delete *CT_it;
                        CT.erase(CT_it);
                        CT_it--;
                        
                        numberOfSaddlesRemoved++;
                        secondStageRemovalDone = true;
                        thisVertexRemoved = true;
                        
                    }
                    
                }
                
                
                if ( (!thisVertexRemoved) && ( (*CT_it)->fromVertices.size() == 1) && ( (*CT_it)->toVertex->index != (*CT_it)->index  ) ) {
                    
                    // Vertex collapse
                    
                    tmp_CT_fromVertex 	= (*CT_it)->fromVertices.front();
                    tmp_CT_toVertex 	= (*CT_it)->toVertex;
                    
                    tmp_El 				= elements[(*CT_it)->index];
                    
                    do {
                        tmp_El->CT_root = tmp_CT_fromVertex;
                        tmp_El = elements[tmp_El->CT_parent];
                    } while( tmp_El->CT_root != tmp_CT_toVertex);
                    
                    std::vector<CT_Vertex*>::iterator tmp_CT_vertex = tmp_CT_toVertex->fromVertices.begin();
                    while( (*tmp_CT_vertex) != *CT_it )
                        tmp_CT_vertex++;
                    tmp_CT_toVertex->fromVertices.erase(tmp_CT_vertex);
                    
                    tmp_CT_toVertex->fromVertices.push_back(tmp_CT_fromVertex);
                    tmp_CT_fromVertex->area = tmp_CT_fromVertex->area + (*CT_it)->area;
                    tmp_CT_fromVertex->toVertex = tmp_CT_toVertex;
                    
                    delete *CT_it;
                    CT.erase(CT_it);
                    CT_it--;
                    
                    numberOfSaddlesRemoved++;
                    secondStageRemovalDone = true;
                    
                }
                
            }
            
        }
        
        contIteration = firstStageRemovalDone || secondStageRemovalDone;
        
    }
    
}


void ContourTree::print_Short() {
    std::cout << "Merge tree size: \t \t " 		<< MT.size() 	<< std::endl;
    std::cout << "Contour tree size: \t \t " 		<< CT.size() 	<< std::endl;
}


uint* ContourTree::export_CT() {
    
    uint* out;
    out = new uint[2*CT.size()]; // a series of edges
    
    for (uint i=0; i<CT.size(); ++i) {
        out[i] 				= CT.at(i)->index;
        out[i+CT.size()] 	= CT.at(i)->toVertex->index;
    }
    
    return out;
}

uint* ContourTree::export_CT_img() {
    
    std::vector<CT_Vertex*>::iterator it;
    Element* tmp_el;
    
    uint N = 1;
    for (uint i=0; i<dimension; ++i)
        N = N*size[i];
    
    uint* out;
    out 	= new uint[N];
    
    for (it = CT.begin(); it != CT.end(); ++it) {
        
        tmp_el = elements[(*it)->index];
        out[(*it)->index] = (*it)->index;
        
        while (tmp_el->CT_parent != (*it)->toVertex->index ) {
            out[tmp_el->CT_parent] = (*it)->index;
            tmp_el	= elements[tmp_el->CT_parent];
        }
        
    }
    
    return out;
    
}

void ContourTree::print_CT(char * ctFileName) {
    std::cout << std::endl;
    std::vector<CT_Vertex*>::iterator it;
    
    std::ofstream out;
    out.open (ctFileName);
    
    for (it = CT.begin(); it != CT.end(); ++it) {
        out		<< (*it)->index << " \t";
        out 		<< (*it)->toVertex->index << " \t";
        out 		<< img[(*it)->index] << " \t";
        out 		<< img[(*it)->toVertex->index] << " \t";
        out		<< (*it)->area;
        out 		<< std::endl;
    }
    
    out.close();
}


void ContourTree::print_CT_img(char * ctImgFileName) {
    
    std::vector<CT_Vertex*>::iterator it;
    Element* tmp_el;
    
    uint N = 1;
    for (uint i=0; i<dimension; ++i)
        N = N*size[i];
    
    double* CT_img;
    CT_img 	= new double[N];
    
    for (it = CT.begin(); it != CT.end(); ++it) {
        
        tmp_el = elements[(*it)->index];
        CT_img[(*it)->index] = img[(*it)->index];
        
        while (tmp_el->CT_parent != (*it)->toVertex->index ) {
            CT_img[tmp_el->CT_parent] = img[(*it)->index];
            tmp_el	= elements[tmp_el->CT_parent];
        }
        
    }
    
    FILE * outImageFile;
    outImageFile = fopen (ctImgFileName, "w");
    
    fwrite (CT_img , sizeof(double), N, outImageFile);
    fclose(outImageFile);
    
    delete[] CT_img;
    
}
