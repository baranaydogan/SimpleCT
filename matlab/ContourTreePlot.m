% ------------------------------------------------------------
% Dogu Baran Aydogan - baran.aydogan@gmail.com
% 21.07.2014
% ------------------------------------------------------------
% 
% "ContourTreePlot.m" plots the contour tree computed by the ContourTree
% function.
%
% - The first number in a vertex box is the index of the pixel/voxel.
% - The second number which is in paranthesis is the intensity (written
% only if ct_img is also provided)
% - The numbers on the edges indicate the number of pixels/voxels on this
% edge
%
% 
% EXAMPLES:
%     [ct ct_img] = ContourTree(img);
%     ContourTreePlot(img, ct);
%     ContourTreePlot(img, ct, ct_img);
%     
% _____________________________________________________________
% WARNING:
% The author does not accept any responsibility or liability for the 
% accuracy of the output since it has not been tested exhaustively. 
% Please use at your own risk.
% 
% _____________________________________________________________
% "ContourTreePlot.m" is a part of SimpleCT.
% 
% SimpleCT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SimpleCT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% _____________________________________________________________


function ContourTreePlot(img, ct, ct_img)

[C,~,ic]                  = unique(ct);
newVertices               = 1:length(C);

tree                      = reshape(newVertices(ic), size(ct));
tree                      = [img(ct(:,1)) tree zeros(size(tree,1),1)];

for i=1:size(tree,1)
    if ( img(ct(i,1)) > img(ct(i,2)) )
        tree(i,4) = 1;
    else
        tree(i,4) = -1;
    end
end

if nargin == 3
    tree = [tree zeros(size(tree,1),1)];
    for i=1:size(tree,1)
        tree(i,5) = nnz(ct_img==ct(i));
    end
end

nodeCount                 = max(max(tree(:,2:3)));
adjacencyMatrix           = zeros( nodeCount );

for i = 1 : nodeCount
    
    fromNode      = i;
    checkToNode   = find(tree(:,2) == i , 1);
    
    if ~isempty(checkToNode)
        
        toNode                  = tree(checkToNode, 3);
        
        direction               = tree(checkToNode, 4);
        
        if size(tree,2) == 5
            weight              = tree(checkToNode, 5);
        end
        
        dispIndex = ct(checkToNode);
        if isempty(dispIndex)
            showInd = '';
        else
            showInd = num2str(dispIndex(1));
        end
        
        if direction == 1
            if size(tree,2) == 5
                adjacencyMatrix(fromNode, toNode) = weight;
                treeIDs{fromNode,1}               = [showInd , ' (', num2str(tree(checkToNode,1)), ')'];
            else
                adjacencyMatrix(fromNode, toNode) = 1;
                treeIDs{fromNode,1}               = showInd;
            end
        elseif direction == -1
            if size(tree,2) == 5
                adjacencyMatrix(toNode, fromNode) = weight;
                treeIDs{fromNode,1}               = [showInd , ' (', num2str(tree(checkToNode,1)), ')'];
            else
                adjacencyMatrix(toNode, fromNode) = 1;
                treeIDs{fromNode,1}               = showInd;
            end
        end
        
    else
        treeIDs{fromNode,1} = num2str(fromNode);
    end
    
    adjacencyMatrix(i,i) = 0;
    
end

if size(tree,2) == 5
    treeGraph = biograph(adjacencyMatrix, treeIDs, 'ShowArrows', 'off', 'ShowWeights', 'on');
else
    treeGraph = biograph(adjacencyMatrix, treeIDs, 'ShowArrows', 'off', 'ShowWeights', 'off');
end

view(treeGraph);
