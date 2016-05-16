function out=beam2(mode,b,c,d,e)
  
% BEAM3 does as listed below. It is an Euler-Bernoulli
% beam/rod/torsion model. 
% Beam properties (bprops) are in the order
% bprops=[E G rho A1 A2 A3 J1 J2 J3 Ixx1 Ixx2 Ixx3 Iyy1 Iyy2 Iyy3]
% For a linear interpolation they are
% bprops=[E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2]
% Note that the linear interpolation a user shortcut
% and results in no less computational effort.
% Third node is in the middle.
% Fourth "node" defines the beam y plane and is actually from the
% points array.
%
%
% Defining beam element properties in wfem input file:
% element properties
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 
%   E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2 
%   E G rho A J Izz Iyy 
%   E G rho A J Izz Iyy sx2 sy2 sz2 srx2 sry2 srz2 distype 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...  
%       mx2 my2 mz2 mrx2 mry2 mrz2 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...  
%       sx2 sy2 sz2 srx2 sry2 srz2 distype 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...  
%       mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 distype 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
%       mx2 my2 mz2 mrx2 mry2 mrz3 mx3 my3 mz3 mrx3 mry3 mrz3 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
%       mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...  
%       mx3 my3 mz3 mrx3 mry3 mrz3 sx3 sy3 sz3 srx3 sry3 srz3 distype
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
%       mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...  
%       lx2 ly2 lz2 lrx2 lry2 lrz2 
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 ...
%       mx2 my2 mz2 mrx2 mry2 mrz2 sx2 sy2 sz2 srx2 sry2 srz2 ...  
%       lx2 ly2 lz2 lrx2 lry2 lrz2 ...  
%       mx3 my3 mz3 mrx3 mry3 mrz3 sx3 sy3 sz3 srx3 sry3 srz3 ...
%       lx3 ly3 lz3 lrx3 lry3 lrz3 
% If the second format is used, a linear interpolation of the
% properties is presumed. If the third format is used, constant
% properties are assumed. In subsequent lines, initial deflections
% can be prescribed, deterministically, or stochastically
% (random). The character \command{m} stands for \emph{mean} value,
% and \command{r} stands for \emph{rotation} value. If a stochastic
% form is used, a distribution type must be prescribes. A
% \command{distype} of \command{0} means normal distribution with
% standard deviation of \command{s} values as prescribe. A
% \command{distype} of \command{-1} means uniform distribution with
% bounds relative to the mean set by the \command{s} values. A
% truncated Gaussian distribution is demonstrated in the last
% line. There the \command{l} values are limits relative to the
% mean (just as the \command{s} values for a uniform
% distribution). Note that the dots illustrate continuation of the
% same line. Continuations of a line using the M\textsc{atlab}
% notation \ldots can only be used in defining element
% properties. Torsional rigidity, $J$, must be less than or equal
% to $Iyy+Izz$ at any given cross section.  
%
% Defining beam3 element in wfem input file:
%   node1 node2 node3 pointnumber materialnumber 
%   node1 node2 pointnumber materialnumber 
%   node1 node2 materialnumber
%
% node2 is automatically created if it is missing.
%
% See wfem.m for more explanation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables (global):
% -------------------
% K       :    Global stiffness matrix
% Ks      :    Global stiffness buckling matrix
% M       :    Global mass matrix
% nodes   :    [x y z] nodal locations
global ismatnewer
global K
global Ks
global M
global nodes % Node locations
global elprops
global element
global points
global Fepsn % Initial strain "forces". 
global lines
global restart
global reload
global curlineno
global DoverL
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NUMOFNODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(mode,'numofnodes')
    out=2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% GENERATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(mode,'generate')
  elnum=c;b;
    element(elnum).nodes=[b(1) b(2)];%n3 recall node 3 is in the
      element(elnum).properties=b(3);
      element(elnum).point=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ISTRAINFORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Here we figure out what the beam properties mean. If you need
% them in a mode, that mode should be in the if on the next line.
if strcmp(mode,'make')|strcmp(mode,'istrainforces')
  elnum=b;element(elnum);
  bnodes=[element(elnum).nodes element(elnum).point];% The point is
                                                     % referred to
                                                     % as node 4
                                                     % below,
                                                     % although it
                                                     % actually
                                                     % calls the
                                                     % array points
                                                     % to get it's location.
  bprops=elprops(element(elnum).properties).a;% element(elnum).properties 
                                              % stores the
                                              % properties number
                                              % of the current
                                              % elnum. elprops
                                              % contains this data.
  
  E=bprops(1);
  G=bprops(2);
  
  rho=bprops(3);

 if length(bprops)==11
   A1=bprops(4);
   A2=bprops(5);
   J1=bprops(6);
   J2=bprops(7);
   Izz1=bprops(8);
   Izz2=bprops(9);
   Iyy1=bprops(10);
   Iyy2=bprops(11);  
  else 
      warndlg('wrong number of material inputs, idiot!')
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% MAKE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beam properties (bprops) are in the order
% bprops=[E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3]
% For a linear beam they are
% bprops=[E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2]

if strcmp(mode,'make')
disp('make')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Define beam node locations for easy later referencing
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  x1=nodes(bnodes(1),1);
  y1=nodes(bnodes(1),2);
  z1=nodes(bnodes(1),3);
  x2=nodes(bnodes(2),1);
  y2=nodes(bnodes(2),2);
  z2=nodes(bnodes(2),3);
  x4=points(bnodes(3),1);
  y4=points(bnodes(3),2);
  z4=points(bnodes(3),3);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Shape functions for higher order beam. 
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Shape functions in matrix polynomial form (polyval style) for bending
  bn1=[2 -3 0 1];
  bn1d=[6 -6 0];
  bn1dd=[12 -6];
  bn2=[1 -2 1 0];
  bn2d=[3 -4 1];
  bn2dd=[6 -4];
  bn3=[-2 3 0 0];
  bn3d=[-6 6 0];
  bn3dd=[-12 6];
  bn4=[1 -1 0 0];
  bn4d=[3 -2 0];
  bn4dd=[6 -2];
  
  % Shape functions in matrix polynomial form (polyval style) for 
  % torsion/rod
  rn1=[-1 1];
  rn1d=[-1];
  rn2=[1 0];
  rn2d=[1];
  numbeamgauss=1; % Number of Gauss points for integration of beam element
  [bgpts,bgpw]=gauss(numbeamgauss)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SOMEONE SCREWED UP %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if numbeamgauss==1,
      bgpts=0
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  kb1=zeros(4,4);
  kb2=kb1;
  l=norm([x2 y2 z2]-[x1 y1 z1]);
  propertynum=num2str(element(elnum).properties);
  % Allowable aspect ratio. I recommend D/l=.1
%   if isempty(DoverL)==1
%     DoverL=.1;
%   end
%   
%   if sqrt(A1*4/pi)/l>DoverL|sqrt(A2*4/pi)/l>DoverL
%     warndlg({['Dimensions of element ' num2str(elnum) ' using properties '...
% 	      propertynum ' are more suitable for a Timoshenko beam.'];...
% 	     'radius divided by length is too large'},...
% 	    'Improper application of element.','replace')
%   end
%   if (Izz1+Iyy1)<(1/2.1*A1^2/pi)|(Izz2+Iyy2)<(1/2.1*A2^2/pi)
%     %2.0 would be exact for a circle
%     warndlg({['Iyy+Izz for properties number' propertynum ' can''t be as '...
% 	      'low as have been given.'];...
% 	     'Nonphysical properties.'},['Impossible cross sectional' ...
% 		    ' properties'],'replace')
%   end
%   slenderness=min([sqrt((Izz1+Iyy1)/A1) sqrt((Izz2+Iyy2)/A2)])/l;
%   if slenderness<.002
%     disp([num2str(elnum) ['is a rediculously thin element. Please' ...
% 		    ' check numbers., and spelling']])
%   end
%   
  Jac=l/2;% Beam Jacobian. valid only if node three is in the
          % middle of the beam
          % Local Bending in x-y plane
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1dd,bgpts(i))/Jac^2;
             polyval(bn2dd,bgpts(i))/Jac;
             polyval(bn3dd,bgpts(i))/Jac^2;
             polyval(bn4dd,bgpts(i))/Jac];
    Izz=polyval(rn1*Izz1+rn2*Izz2,bgpts(i));%these should
                                                     %be called Izz
    kb1=kb1+bgpw(i)*beamsfs*beamsfs'*Izz*E*Jac
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
  numrodgauss=3;
  [rgpts,rgpw]=gauss(numrodgauss);
  krod=zeros(2,2);
  ktor=zeros(2,2);
  for i=1:numrodgauss
    rodsfs=[polyval(rn1d,rgpts(i))/Jac;
            polyval(rn2d,rgpts(i))/Jac];
            
    if (J1>(Iyy1+Izz1))|(J2>(Iyy2+Izz2))
      if (J1>(Iyy1+Izz1))
	disp('WARNING: J1 must be <= Iyy1+Izz1')
      end
      if (J2>(Iyy2+Izz2))
	disp('WARNING: J2 must be <= Iyy2+Izz2')
      
      disp(['Error in element properties number '... 
	    num2str(element(elnum).properties) ...
	    'used by element ' num2str(elnum) ' on line'...
	    num2str(element(elnum).lineno) '.'])
    end
    J=polyval(rn1*J1+rn2*J2,bgpts(i));
    A=polyval(rn1*A1+rn2*A2,bgpts(i));
    krod=krod+rgpw(i)*rodsfs*rodsfs'*A*E*Jac;
    ktor=ktor+rgpw(i)*rodsfs*rodsfs'*J*G*Jac;
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % Derivation of Mass matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  numbeamgauss=numbeamgauss+3;
  [bgpts,bgpw]=gauss(numbeamgauss);
  mb1=zeros(4,4);
  % Local Bending in x-y plane
  for i=1:numbeamgauss
    beamsfs=[polyval(bn1,bgpts(i));
             polyval(bn2,bgpts(i))*Jac;
             polyval(bn3,bgpts(i));
             polyval(bn4,bgpts(i))*Jac];
    A=polyval(rn1*A1+rn2*A2,bgpts(i));
    mb1=mb1+bgpw(i)*beamsfs*beamsfs'*rho*A*Jac;%pause
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
  numrodgauss=numrodgauss+1;
  [rgpts,rgpw]=gauss(numrodgauss);
  mrod=zeros(2,2);
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
  % stiffness matrix
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

  R1=([x2 y2 z2]-[x1 y1 z1]);
  lam1=R1/norm(R1);
  R2=([x4 y4 z4]-[x1 y1 z1]);
  R2perp=R2-dot(R2,lam1)*lam1;
  udirec=0;
  while norm(R2perp)<10*eps
    udirec=udirec+1;
    %disp('oops')
    %pause
    [minval,minloc]=min(lam1);
    R2perp=zeros(1,3);
    R2perp(udirec)=1;
    R2perp=R2perp-dot(R2perp,lam1)*lam1;
  end
  lam2=R2perp/norm(R2perp);
  lam3=cross(lam1,lam2);
  lamloc=[lam1;lam2;lam3];
  lam=sparse(12,12);
  lam(1:3,1:3)=lamloc;
  lam(4:6,4:6)=lamloc;
  lam(7:9,7:9)=lamloc;
  lam(10:12,10:12)=lamloc;
%   lam(13:15,13:15)=lamloc;
%   lam(16:18,16:18)=lamloc;
  
% $$$     lam=[lamloc z z z z z;
% $$$          z lamloc z z z z;
% $$$          z z lamloc z z z;
% $$$          z z z lamloc z z;
% $$$          z z z z lamloc z;
% $$$          z z z z z lamloc];
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
kg

  K(indices,indices)=K(indices,indices)+kg;
  M(indices,indices)=M(indices,indices)+mg;

  % At this point we also know how to draw the element (what lines
  % and surfaces exist). For the beam3 element, 2 lines are
  % appropriate.
  numlines=size(lines,1);
  lines(numlines+1,:)=[bn1 bn2];
  %lines(numlines+2,:)=[bn3 bn2];
  
%diag(M)
% elseif strcmp(mode,'istrainforces')
%   %disp('beam3 initial strain forces')
%   % We need to have the stiffness matrix and the coordinate roation matrix.
%   
%   k=element(elnum).k;
%   lam=element(elnum).lambda;
%   kg=lam'*k*lam;
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %
%   % Add initial deflections
%   %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   if exist('distype')==0
%     distype=0;
%   end
%   %disp('before loop')
%   %bprops
%   %length(bprops)
%   if length(bprops)>15|length(bprops)==14
%     %disp('in loop')
%     if distype==0
%       %disp('need to check for limits, then apply them/make them')
%       %normal distribution
%       if ~exist('ldx2')
%         ldx2=inf;
%         ldy2=inf;
%         ldz2=inf;
%         ldrx2=inf;
%         ldry2=inf;
%         ldrz2=inf;
%         ldx3=inf;
%         ldy3=inf;
%         ldz3=inf;
%         ldrx3=inf;
%         ldry3=inf;
%         ldrz3=inf;
%       end
%       dev=randn(1)*sdx2;
%       while dev>ldx2
%         dev=randn(1)*sdx2;
%       end
%       dx2=mdx2+dev;
%       dev=randn(1)*sdy2;
%       while dev>ldy2
%         dev=randn(1)*sdy2;
%       end
%       dy2=mdy2+dev;
%       dev=randn(1)*sdz2;
%       while dev>ldz2
%         dev=randn(1)*sdz2;
%       end
%       dz2=mdz2+dev;
%       dev=randn(1)*sdrx2;
%       while dev>ldrx2
%         dev=randn(1)*sdrx2;
%       end
%       drx2=mdrx2+dev;
%       dev=randn(1)*sdry2;
%       while dev>ldry2
%         dev=randn(1)*sdry2;
%       end
%       dry2=mdry2+dev;
%       dev=randn(1)*sdrz2;
%       while dev>ldrz2
%         dev=randn(1)*sdrz2;
%       end
%       drz2=mdrz2+dev;
%     elseif abs(distype)==1
%       %uniform distribution
%       %disp('uniform distribution')
%       dx2=mdx2+(rand(1)-.5)*sdx2*2;%pause
%       dy2=mdy2+(rand(1)-.5)*sdy2*2;
%       dz2=mdz2+(rand(1)-.5)*sdz2*2;
%       drx2=mdrx2+(rand(1)-.5)*sdrx2*2;
%       dry2=mdry2+(rand(1)-.5)*sdry2*2;
%       drz2=mdrz2+(rand(1)-.5)*sdrz2*2;
%     else
%       disp(['Invalid distribution type ' num2str(distype) '.'])
%     end
%     if exist('mdx3')==0
%       %Calculate 'smooth' displacements of node 3 if they are
%       %undefined
% 
%       node3disps=-k(13:18,13:18)\k(13:18,7:12)*[dx2;dy2;dz2;drx2;dry2;drz2];
%       dx3=node3disps(1);
%       dy3=node3disps(2);
%       dz3=node3disps(3);
%       drx3=node3disps(4);
%       dry3=node3disps(5);
%       drz3=node3disps(6);
%     else
%       if distype==0
%         %normal distribution
%         dev=randn(1)*sdx3;
%         while dev>ldx3
%           dev=randn(1)*sdx3;
%         end
%         dx3=mdx3+dev;
%         dev=randn(1)*sdy3;
%         while dev>ldy3
%           dev=randn(1)*sdy3;
%         end
%         dy3=mdy3+dev;
%         dev=randn(1)*sdz3;
%         while dev>ldz3
%           dev=randn(1)*sdz3;
%         end
%         dz3=mdz3+dev
%         dev=randn(1)*sdrx3;
%         while dev>ldrx3
%           dev=randn(1)*sdrx3;
%         end
%         drx3=mdrx3+dev;
%         dev=randn(1)*sdry3;
%         while dev>ldry3
%           dev=randn(1)*sdry3;
%         end
%         dry3=mdry3+dev;
%         dev=randn(1)*sdrz3;
%         while dev>ldrz3
%           dev=randn(1)*sdrz3;
%         end
%         drz3=mdrz3+dev;
%       elseif abs(distype)==1
%         %uniform distribution
%         dx3=mdx3+(rand(1)-.5)*sdx3*2;
%         dy3=mdy3+(rand(1)-.5)*sdy3*2;
%         dz3=mdz3+(rand(1)-.5)*sdz3*2;
%         drx3=mdrx3+(rand(1)-.5)*sdrx3*2;
%         dry3=mdry3+(rand(1)-.5)*sdry3*2;
%         drz3=mdrz3+(rand(1)-.5)*sdrz3*2;
%       else
%         disp(['Invalid distribution type ' num2str(distype) '.'] )    
%       end    
%     end
%     %disp('hello'),pause
%     Fistrain=k*[0;0;0;0;0;0;dx2;dy2;dz2;drx2;dry2;drz2;dx3;dy3;dz3;drx3;dry3;drz3];
%     % Need to transform Fistrain into global coordinates
%     Fistrain=lam'*Fistrain;
%     %  mg=lam'*m*lam;
%   
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %
%     % Assembling matrices into global matrices
%     %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     bn1=bnodes(1);bn2=bnodes(2);bn3=bnodes(3);
%     indices=[bn1*6+(-5:0) bn2*6+(-5:0) bn3*6+(-5:0)] ;
%     
%     if length(Fepsn)<max(indices);
%       Fepsn(max(indices))=0;
%     end
%     Fepsn(indices)=Fepsn(indices)+Fistrain;
%   end
%   
% elseif strcmp(mode,'draw')
% elseif strcmp(mode,'buckle')
  end
end

