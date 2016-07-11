function out=BRICK_Ado(mode,b,c,d,e)
 
global K
global M
global nodes
global elprops
global element
global points
global lines
global surfs

%
% Variables (local):
% ------------------
% bnodes  :    node/point numbers for actual beam nodes 1-2-3 and point
% k       :    stiffness matrix in local coordiates
% kg      :    stiffness matrix rotated into global coordinates
% m       :    mass matrix in local coordiates
% mg      :    mass matrix rotated into global coordinates
% 

out=0;
if strcmp(mode,'numofnodes')
    % This allows a code to find out how many nodes this element has
    out=8;
end
if strcmp(mode,'generate')
  elnum=c;%When this mode is called, the element number is the 3rd
          %argument. 
          %The second argument (b) is the element
          %definition. For this element b is
          %node1 node2 point(for rotation) and material#
          
  if length(b)==9
      element(elnum).nodes=b(1:8);
      element(elnum).properties=b(9);
      %element(elnum).point=b(3);
  else 
	  b
      %There have to be four numbers on a line defining the
      %element. 
      warndlg(['Element ' num2str(elnum) ' on line ' ...
               num2str(element(elnum).lineno) ' entered incorrectly.'], ...
              ['Malformed Element'],'modal')
      return
  end
 
end

% Here we figure out what the beam properties mean. If you need
% them in a mode, that mode should be in the if on the next line.

if strcmp(mode,'make')||strcmp(mode,'istrainforces')
  elnum=b;% When this mode is called, the element number is given
          % as the second input.
          
  bnodes=[element(elnum).nodes];% The point is referred to
  % as node 4 below, although it actually calls the array points to get its
  % location. Its not really a node, but just a point that helps define
  % orientation. Your element may not need such a reference point.
  
  bprops=elprops(element(elnum).properties).a;% element(elnum).properties 
  % stores the properties number of the current elnum. elprops contains 
  % this data. This is precisely the material properties line in an
  % array. You can pull out any value you need for your use. 
  
  if length(bprops)==3
      % Should only need 3 properties to generate the brick element
      E=bprops(1); mu=bprops(2); rho=bprops(3);
      
  else
      warndlg(['The number of material properties set for ' ...
               'this element (' num2str(length(bprops)) ') isn''t ' ...
               'appropriate for a brick element. '      ...
               'Please refer to the manual.'],...
              'Bad element property definition.','modal');
  end
end

if strcmp(mode,'make') % Define BRICK node locations for easy later referencing
  
  x1=nodes(bnodes(1),1);  % Location of node 1 x1 y1 z1
  y1=nodes(bnodes(1),2);
  z1=nodes(bnodes(1),3);
  x2=nodes(bnodes(2),1);  % Location of node 2 x2 y2 z2
  y2=nodes(bnodes(2),2);
  z2=nodes(bnodes(2),3);
  x3=nodes(bnodes(3),1);  % Location of node 3 x3 y3 z3
  y3=nodes(bnodes(3),2);
  z3=nodes(bnodes(3),3);
  x4=nodes(bnodes(4),1);  % Location of node 4 x4 y4 z4
  y4=nodes(bnodes(4),2);
  z4=nodes(bnodes(4),3);
  x5=nodes(bnodes(5),1);  % Location of node 5 x5 y5 z5
  y5=nodes(bnodes(5),2);
  z5=nodes(bnodes(5),3); 
  x6=nodes(bnodes(6),1);  % Location of node 6 x6 y6 z6
  y6=nodes(bnodes(6),2);
  z6=nodes(bnodes(6),3); 
  x7=nodes(bnodes(7),1);  % Location of node 7 x7 y7 z7
  y7=nodes(bnodes(7),2);
  z7=nodes(bnodes(7),3); 
  x8=nodes(bnodes(8),1);  % Location of node 8 x8 y8 z8
  y8=nodes(bnodes(8),2);
  z8=nodes(bnodes(8),3); 
  %
  
  global_nodes = [[x1 y1 z1]; [x2 y2 z2]; [x3 y3 z3]; [x4 y4 z4];...
      [x5 y5 z5]; [x6 y6 z6]; [x7 y7 z7]; [x8 y8 z8]];

  % Number of Gauss points for integration of BRICK element
  % numbeamgauss=2; [bgpts,bgpw]=gauss(numbeamgauss);
  
  num_gauss=2
  [int_p, int_w] = gauss(num_gauss)
  intPts=zeros(num_gauss^3,3);
  intWts=zeros(num_gauss^3,3);
  index=0;
  
  for i=1:num_gauss
      for j=1:num_gauss
          for k=1:num_gauss
              index=index+1
              intPts(index,:) = [int_p(i) int_p(j) int_p(k)];
              intWts(index,:) = [int_w(i) int_w(j) int_w(k)];
          end
      end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % Derivation of Stiffness and Mass matrix
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  B=zeros(6,24);
  Ke=zeros(24,24);
  Me=zeros(24,24);
  id=[1 4 7 10 13 16 19 22];
  
  Em=getE(E,mu);
  
  % Loop to construct BRICK stiffness matrix
  for p=1:num_gauss^3
      r = intPts(p,1);
      s = intPts(p,2);
      t = intPts(p,3);
      
      Ne=getN(r,s,t);
      dN=getdN(r,s,t);
      
      J=dN*global_nodes;
      Jinv=J\eye(3);
      JDet = det(J);
      
      for q=1:8
          dN_i = dN(:,q);
          
          Bi=[Jinv(1,:)*dN_i 0 0;
              0 Jinv(2,:)*dN_i 0;
              0 0 Jinv(3,:)*dN_i;
              Jinv(2,:)*dN_i Jinv(1,:)*dN_i 0;
              0 Jinv(3,:)*dN_i Jinv(2,:)*dN_i;
              Jinv(3,:)*dN_i 0 Jinv(1,:)*dN_i];
          
          B(1:end, id(q):id(q)+2) = Bi(1:end, 1:end);
      end
      
      Ki = prod(intWts(p,1:end))*JDet*(B'*Em*B);
      Ke=Ke+Ki;
      
      Mi = rho*prod(intWts(p,1:end))*JDet*(Ne'*Ne);
      Me=Me+Mi;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  element(elnum).m=Me;
  element(elnum).k=Ke;

  indices = zeros(1,24);
  for w = 1:8
      currentNode = element(elnum).nodes(w);
      
      indices(3*w-2:3*w) = 3*currentNode-2:3*currentNode;
      
  end

  K(indices,indices)=K(indices,indices)+Ke;
  M(indices,indices)=M(indices,indices)+Me;

%   % Connecting node information
%   numlines=size(lines,1);
%   lines(numlines+1,:)=[bn1 bn2];
  
  end
end
