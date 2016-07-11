function out=beam2exampleMCk(mode,b,c,d,e)
  
% BEAM2 does as listed below. It is an Euler-Bernoulli
% beam/rod/torsion model. 
% Beam properties (bprops) are in the order
% bprops=[E G rho A1 A2 J1 J2 Ixx1 Ixx2 Iyy1 Iyy2]
% Third node is in the middle.
% Fourth "node" defines the beam y plane and is actually from the
% points array.
%
% Defining beam element properties in wfem input file:
% element properties
%   E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2 
% Torsional rigidity, $J$, must be less than or equal
% to $Iyy+Izz$ at any given cross section.  
%
% Defining beam2 element in wfem input file:
%   node1 node2 node3 pointnumber materialnumber 

%
% See wfem.m for more explanation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables (global):
% -------------------
% K       :    Global stiffness matrix
% M       :    Global mass matrix
% nodes   :    [x y z] nodal locations
global K
global M
global nodes
global elprops
global element
global points
global lines
%global DoverL

%
% Variables (local):
% ------------------
% bnodes  :    node/point numbers for actual beam nodes 1-2-3 and point
% k       :    stiffness matrix in local coordiates
% kg      :    stiffness matrix rotated into global coordinates
% m       :    mass matrix in local coordiates
% mg      :    mass matrix rotated into global coordinates
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Joseph C. Slater, 7/26/2002.
% joseph.slater@wright.edu
out=0;
if strcmp(mode,'numofnodes')
    % This allows a code to find out how many nodes this element has
    out=2;
end
if strcmp(mode,'generate')
  elnum=c;%When this mode is called, the element number is the 3rd
          %argument. 
          %The second argument (b) is the element
          %definition. For this element b is
          %node1 node2 point(for rotation) and material#
          
  if length(b)==4
      element(elnum).nodes=b(1:2);
      element(elnum).properties=b(4);
      element(elnum).point=b(3);
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
          
  bnodes=[element(elnum).nodes element(elnum).point];% The point is referred to
  % as node 4 below, although it actually calls the array points to get its
  % location. Its not really a node, but just a point that helps define
  % orientation. Your element may not need such a reference point.
  
  bprops=elprops(element(elnum).properties).a;% element(elnum).properties 
  % stores the properties number of the current elnum. elprops contains 
  % this data. This is precisely the material properties line in an
  % array. You can pull out any value you need for your use. 
  
  if length(bprops)==11
      % Beam properties (bprops) are in the order
      % bprops=[E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2]
      E=bprops(1); G=bprops(2); rho=bprops(3);
      A1=bprops(4); A2=bprops(5);
      J1=bprops(6); J2=bprops(7);
      Izz1=bprops(8); Izz2=bprops(9); Iyy1=bprops(10); Iyy2=bprops(11);

  else
      warndlg(['The number of material properties set for ' ...
               'this element (' num2str(length(bprops)) ') isn''t ' ...
               'appropriate for a beam2 element. '      ...
               'Please refer to the manual.'],...
              'Bad element property definition.','modal');
  end
end

if strcmp(mode,'make') % Define beam node locations for easy later referencing
  
  x1=nodes(bnodes(1),1);  % Location of node 1 x1 y1 z1
  y1=nodes(bnodes(1),2);
  z1=nodes(bnodes(1),3);
  x2=nodes(bnodes(2),1);  % Location of node 2 x2 y2 z2
  y2=nodes(bnodes(2),2);
  z2=nodes(bnodes(2),3);

  %
  % Shape functions for higher order beam. 
  % Shape functions in matrix polynomial form (polyval style) for bending
  bn1 =  1/4*[ 1 0 -3 2]; 
  % bn1d =  1/4*[ 3 0 -3 ]; 
  bn1dd =  [ 3/2 0 ]; 
  bn2 =  1/4*[ 1 -1 -1 1 ]; 
  %bn2d =  1/4*[ 3 -2 -1 ];  
  bn2dd =  [ 3/2 -1/2 ]; 
  bn3 =  1/4*[ -1 0 3 2 ];  
  %bn3d = 1/4*[ -3 0 3 ];  
  bn3dd = [ -3/2 0 ];     
  bn4 =  1/4*[ 1 1 -1 -1 ];  
  %bn4d = 1/4*[ 3 2 -1 ]; 
  bn4dd = [ 3/2 1/2 ]; 
  
  % Shape functions in matrix polynomial form (polyval style) for 
  % torsion/rod
  rn1= [-0.5 0.5];
  rn1d= -0.5;
  rn2= [0.5 0.5];
  rn2d= 0.5;

  % Number of Gauss points for integration of beam element
  numbeamgauss=5; [bgpts,bgpw]=gauss(numbeamgauss);
  
  % For this beam, 2 nodes, 2DOF each, is a 4 by 4 matrix. 
  kb1=zeros(4,4); kb2=kb1;

  propertynum=num2str(element(elnum).properties);
  
  Jac=l/2;% Beam Jacobian. 
  % Local Bending in x-y plane
  for i=1:numbeamgauss
    % evaluating second derivatives of shape functions to use in
    % generating stiffness matrix. (at gauss point)
    beamsfs=[polyval(bn1dd,bgpts(i))/Jac^2;
             polyval(bn2dd,bgpts(i))/Jac;
             polyval(bn3dd,bgpts(i))/Jac^2;
             polyval(bn4dd,bgpts(i))/Jac];
    %Find Izz at Gauss point
    Izz=polyval(rn1*Izz1+rn2*Izz2,bgpts(i));
    %This is the Gauss integration part. 
    kb1=kb1+bgpw(i)*beamsfs*beamsfs'*Izz*E*Jac;
  end
  % Local Bending in x-z plane
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1dd,bgpts(i))/Jac^2;
             -polyval(bn2dd,bgpts(i))/Jac;
             polyval(bn3dd,bgpts(i))/Jac^2;
             -polyval(bn4dd,bgpts(i))/Jac];
    Iyy=polyval(rn1*Iyy1+rn2*Iyy2,bgpts(i));
    kb2=kb2+bgpw(i)*beamsfs*beamsfs'*Iyy*E*Jac;
  end
  
  % Local Extension in x, torsion about x
  numrodgauss=3;% Number of points to use for gauss point integration
  [rgpts,rgpw]=gauss(numrodgauss);
  krod=zeros(2,2);
  ktor=zeros(2,2);
  for i=1:numrodgauss
    rodsfs=[polyval(rn1d,rgpts(i))/Jac;
            polyval(rn2d,rgpts(i))/Jac];

    J=polyval(rn1*J1+rn2*J2,bgpts(i));% J at gauss point.
    A=polyval(rn1*A1+rn2*A2,bgpts(i));% A at gauss point
    
    % Since the shape functions and Gauss points are the same, we are doing
    % the rod and torsion rod together. 
    krod=krod+rgpw(i)*rodsfs*rodsfs'*A*E*Jac;
    ktor=ktor+rgpw(i)*rodsfs*rodsfs'*J*G*Jac;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % Derivation of Mass matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  numbeamgauss=numbeamgauss+3; %Need more gauss points for the mass
                               %matrix. 
  [bgpts,bgpw]=gauss(numbeamgauss);
  mb1=zeros(4,4); %initialize empty mass matrix
  % Local Bending in x-y plane
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1,bgpts(i));
             polyval(bn2,bgpts(i))*Jac;
             polyval(bn3,bgpts(i));
             polyval(bn4,bgpts(i))*Jac];
    A=polyval(rn1*A1+rn2*A2,bgpts(i));
    mb1=mb1+bgpw(i)*beamsfs*beamsfs'*rho*A*Jac;
  end
  
  % Local Bending in x-z plane
  mb2=zeros(4,4);
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1,bgpts(i));
             -polyval(bn2,bgpts(i))*Jac;
             polyval(bn3,bgpts(i));
             -polyval(bn4,bgpts(i))*Jac];
    A=polyval(rn1*A1+rn2*A2,bgpts(i));
    mb2=mb2+bgpw(i)*beamsfs*beamsfs'*rho*A*Jac;
  end
  
  % Local Extension in x, torsion about x
  numrodgauss=numrodgauss+1; %Need more gauss points for the mass
                             %matrix. 
  [rgpts,rgpw]=gauss(numrodgauss);
  mrod=zeros(2,2); %initialize empty mass matrix
  mtor=zeros(2,2);
  for i=1:numrodgauss
    rodsfs=[polyval(rn1,rgpts(i));
            polyval(rn2,rgpts(i))];
    J=polyval(rn1*(Iyy1+Izz1)+rn2*(Iyy2+Izz2),bgpts(i));
    A=polyval(rn1*A1+rn2*A2,bgpts(i));
    mrod=mrod+rgpw(i)*rodsfs*rodsfs'*A*rho*Jac;
    mtor=mtor+rgpw(i)*rodsfs*rodsfs'*J*rho*Jac;
  end
  
  % Assembling each stiffness matrix into the complete elemental 
  % stiffness matrix. We're just telling the sub-elements to be put
  % into the correct spots for the total element. 
  k=zeros(12,12);
  k([2 6 8 12],[2 6 8 12])=kb1;
  k([3 5 9 11],[3 5 9 11])=kb2;
  k([1 7],[1 7])=krod;
  k([4 10],[4 10])=ktor;
  
  % Assembling each mass matrix into the complete elemental 
  % mass matrix
  m=zeros(12,12);
  m([2 6 8 12],[2 6 8 12])=mb1;
  m([3 5 9 11],[3 5 9 11])=mb2;
  m([1 7],[1 7])=mrod;
  m([4 10],[4 10])=mtor;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Coordinate rotations
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  R1=([x2 y2 z2]-[x1 y1 z1]);% Vector along element
  lam1=R1/norm(R1);% Unit direction
  R2=([x4 y4 z4]-[x1 y1 z1]);% Unit direction to point
  R2perp=R2-dot(R2,lam1)*lam1;% Part of R2 perpendicular to lam1
  udirec=0;
  while norm(R2perp)<10*eps
    udirec=udirec+1;
    [minval,minloc]=min(lam1);
    R2perp=zeros(1,3);
    R2perp(udirec)=1;
    R2perp=R2perp-dot(R2perp,lam1)*lam1;
  end
  %Make the unit direction vectors for rotating and put them in the
  %rotation matrix. 
  lam2=R2perp/norm(R2perp);
  lam3=cross(lam1,lam2);
  lamloc=[lam1;lam2;lam3];
  lam=sparse(12,12);
  lam(1:3,1:3)=lamloc;
  lam(4:6,4:6)=lamloc;
  lam(7:9,7:9)=lamloc;
  lam(10:12,10:12)=lamloc;
  
  % Updating global vars
  element(elnum).lambda=lam;
  element(elnum).m=m;
  element(elnum).k=k;

  kg=lam'*k*lam;
  mg=lam'*m*lam;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  bn1=bnodes(1);bn2=bnodes(2);%bn3=bnodes(3);
  indices=[bn1*6+(-5:0) bn2*6+(-5:0)];% bn3*6+(-5:0)] ;


  K(indices,indices)=K(indices,indices)+kg;
  M(indices,indices)=M(indices,indices)+mg;

  % Connecting node information
  numlines=size(lines,1);
  lines(numlines+1,:)=[bn1 bn2];
  
end
