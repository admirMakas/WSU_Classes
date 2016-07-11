clear all
close all
clc
%input maker
%makes a rectangular pyramid, like an incan pyramid or a ziggurat
levels=5;% Levels? 1, 2, 3, 4, any positive integer'
epl=4;%('Elements per level? 1, 4, 9, 16 must be perfect square, elements per level')
h=2; %height
bw=.5; %base width 
tw=bw; %top width
E=200e9; %young's/elasticity modulus
v=.3; %poisson's ratio
G=8000; %shear modulus
force=6.2500e+09; %can only apply one force at the very last node with this
fdof=3; %degree of freedom the force is acting in
%assuming equal sided rectagular pyramid

sidel=epl^.5+1; %nodes per side
total_elements=levels*epl; %total elements
total_nodes=sidel^2*(levels+1); %total nodes


glat=linspace(-1,1,sidel); %lateral grade
ghor=linspace(bw/2,tw/2,levels+1);
gh=linspace(0,h,levels+1) ;

node_make=zeros(sidel,sidel,levels+1);
node_number=node_make;
xco=node_make;
yco=node_make;
zco=node_make;

%check
if numel(node_make)==total_nodes
    disp('Correct number of nodes')
else
    disp('Incorrect number of nodes')
end

j=[1 2 4 3 5 6 8 7];

%% making node coordinates

k=1;
for plane=1:levels+1
    for row=1:sidel
        for column=1:sidel
          node_number(column,row,plane)=k;
          k=k+1;
          zco(column,row,plane)=gh(plane);
          xco(column,row,plane)=glat(column)*ghor(plane);
          yco(column,row,plane)=glat(row)*ghor(plane);
        end
    end
end

node_number;
vnn=node_number(1:end)';
xc=xco(1:end)';
yc=yco(1:end)';
zc=zco(1:end)';
all_nodes=[vnn xc yc zc];

% plot3(xc,yc,zc,'o')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid on
% text(xc,yc,zc,num2str(node_number(1:end)'))


%% picking elements 

k=1;
for plane=1:levels
    for row=1:sidel-1
        for column=1:sidel-1
         % properties=1
          picked=node_number(column:1+column,row:1+row,plane:1+plane); 
          elementp(k,:)=picked(j);
          k=k+1;
          
        end
    end
end

for enum=1:total_elements;
all_elements(enum,:)=[elementp(enum,:) 1];
end

%irrelevant
%element(#).nodes
% node_number

% for elnum=1:total_elements
% front=element(elnum).nodes(1:4);
% back=element(elnum).nodes(5:8);
% top=element(elnum).nodes([3 4 8 7]);
% bottom=element(elnum).nodes([1 5 6 2]);
% left=element(elnum).nodes([2 3 7 6]);
% right=element(elnum).nodes([1 5 8 4]);
% 
% element(elnum).surfacenodes=[front;back;top;bottom;left;right];
% end

% %restrain these nodes
restrain_nodes=node_number(1:numel(node_number(:,:,1)))';

for n=1:length(restrain_nodes)
   all_restraints(n,:)=[restrain_nodes(n) 0 0 1];
end

%% writing to file

mydata{1}='element properties';
mydata{end+1,1}=num2str([E v G]);
mydata{end+1,1}='';
mydata{end+1,1}='BRICK_Ado2 elements';
for k=1:total_elements;
mydata{end+1,1}=num2str(all_elements(k,:));
end
mydata{end+1,1}='';
mydata{end+1,1}='nodes';
for k=1:length(vnn);
mydata{end+1,1}=num2str(all_nodes(k,:));
end
% out of nodes, elements, and properties
%moving on to points, fixation, loads
mydata{end+1,1}='';
mydata{end+1,1}='points';
mydata{end+1,1}='1 1 1 1';
mydata{end+1,1}='';
mydata{end+1,1}='fix surfaceball';
mydata{end+1,1}='1 0 1 0';
mydata{end+1,1}='1 1 0 0';
mydata{end+1,1}='2 0 1 0';
mydata{end+1,1}='2 1 0 0';
for k=1:length(restrain_nodes);
mydata{end+1,1}=num2str(all_restraints(k,:));
end
mydata{end+1,1}='';
%force
mydata{end+1,1}='load';
npl=numel(node_number(:,:,end)); %nodes per level
for k=(total_nodes-npl+1):total_nodes;
mydata{end+1,1}=num2str([node_number(k) fdof force/npl]);
end
%postprocessing stuff
mydata{end+1,1}='';
mydata{end+1,1}='actions';
mydata{end+1,1}='staticanalysis';
mydata{end+1,1}='post_process';
mydata{end+1,1}='end';
mydata{end+1,1}='';
fileID=fopen('pyramidinput.txt','w');
fprintf(fileID,'%s\n',mydata{:});
fclose(fileID);

d_should_be=force*h/(bw^2*E)

