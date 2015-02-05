/**************************************************************************
 *
 * -----------------------------------------------------------------------
 * Dogu Baran Aydogan - baran.aydogan@gmail.com
 * 21.07.2014
 * -----------------------------------------------------------------------
 *
 *************************************************************************
 * "ContourTree" is a wrapper for the SimpleCT library for Matlab
 *
 *
 * EXAMPLES:
 *      ct          = ContourTree(I);
 *      ct          = ContourTree(I, INTENSITY_THRESH, AREA_THRESH);
 *      [ct ct_img] = ContourTree(I);
 *      [ct ct_img] = ContourTree(I, INTENSITY_THRESH, AREA_THRESH);
 *
 *
 * INPUTS:
 *
 *      I                   2D or 3D input image.
 *
 *      INTENSITY_THRESH    Intensity threshold for pruning
 *      AREA_THRESH         Area/volume threshold for pruning
 *
 * OUTPUT:
 *
 *      ct                  A matrix of contour tree edge
 *      ct_img              An image for the contour tree
 *
 * WARNING:
 * The author does not accept any responsibility or liability for the
 * accuracy of the output since it has not been tested exhaustively.
 * Please use at your own risk.
 *
 *
 * ____________________________________________________________________
 * "ContourTree.cpp" is a part of SimpleCT.
 *
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
 *
 *************************************************************************/

#include <mex.h>
#include <vector>
#include <math.h>
#include <omp.h>
#include "SimpleCT.h"

// Compilation commands
//
// For Linux
// mex ContourTree.cpp SimpleCT.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
//
// For Windows (not tested)
// mex  ContourTree.cpp SimpleCT.cpp COMPFLAGS="$COMPFLAGS /openmp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // *********  Abort if the function is not called properly **********//
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("MATLAB:ContourTree:rhs", "Aborting: This function requires at least 1 input. \n\n");
    }
    if (nrhs > 3) {
        mexErrMsgIdAndTxt("MATLAB:ContourTree:rhs", "Aborting: This function requires at most 3 inputs. \n\n");
    }
    
    double *img     = mxGetPr(prhs[0]);                 // Input image
    uint dim_img    = mxGetNumberOfDimensions(prhs[0]); // Dimension of image
    const mwSize *dims_img;                             // Image size
    dims_img        = mxGetDimensions(prhs[0]);
    
    int numberOfElements;
    
    std::vector<uint> size;             // SimpleCT requires the size of image as a vector input
    size.push_back(dims_img[0]);        // These lines convert MATLABS mwSize type to uint vector
    size.push_back(dims_img[1]);
    numberOfElements = dims_img[0]*dims_img[1];
    if (dim_img > 2) {
        size.push_back(dims_img[2]);
        numberOfElements = numberOfElements*dims_img[2];
    }
    
    double intThresh;   //Variables for pruning parameters
    double areaThresh;
    
    if (nrhs == 3) {
        intThresh   = mxGetScalar(prhs[1]);             // Read intensity threshold
        areaThresh  = mxGetScalar(prhs[2]);             // Read area/volume threshold
    } else {
        intThresh   = 0; // Default pruning parameters
        areaThresh  = 0;
    }
    
    
    ContourTree ct(img,size);                   // Compute merge tree
    ct.prune(areaThresh,intThresh);             // Prune and create the contour tree
    
    uint* ct_matrix;
    ct_matrix           = ct.export_CT();       // Export contour tree as a graph
    
    
    int numberOfEdges   = ct.getCT()->size();
    
    plhs[0]             = mxCreateDoubleMatrix(numberOfEdges, 2, mxREAL);
    double*  ct_out     = mxGetPr(plhs[0]);
    for (int i = 0; i < 2*numberOfEdges; ++i)
        ct_out[i]       = ct_matrix[i]+1;
    delete[] ct_matrix;                         // Don't forget to delete, this memory was allocated by SimpleCT
    
    uint*   ct_img;
    double* ct_img_out;
    
    if (nlhs == 2) {
        ct_img              = ct.export_CT_img();   // Export contour tree image
        plhs[1]             = mxCreateNumericArray(dim_img, dims_img, mxDOUBLE_CLASS, mxREAL);
        ct_img_out          = mxGetPr(plhs[1]);
        
        for (int i = 0; i < numberOfElements; ++i)
            ct_img_out[i]       = double(ct_img[i])+1;
        
        delete[] ct_img;                         // Don't forget to delete, this memory was allocated by SimpleCT
    }
    
}
