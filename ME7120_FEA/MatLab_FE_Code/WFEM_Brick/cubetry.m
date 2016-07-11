global displacements
global nodes
global element
global surfaces
global total_elements
global elnum
global X

%outer_limits=max(nodes)
%total_elements=3

figure
scale = 30;
for elnum=1:total_elements
    
    x=zeros(4,6);
    y=zeros(4,6);
    z=zeros(4,6);
    
    dx=displacements(1,:)*scale;
    dy=displacements(2,:)*scale;
    dz=displacements(3,:)*scale;
    
    for j=1:4; %for node
        for i=1:6; %for face
            
            x(j,i)=nodes(element(elnum).surfacenodes(i,j),1);
            y(j,i)=nodes(element(elnum).surfacenodes(i,j),2);
            z(j,i)=nodes(element(elnum).surfacenodes(i,j),3);
            
            xdx(j,i)=x(j,i)+dx(element(elnum).surfacenodes(i,j));
            ydy(j,i)=y(j,i)+dy(element(elnum).surfacenodes(i,j));
            zdz(j,i)=z(j,i)+dz(element(elnum).surfacenodes(i,j));
            
        end
    end

for i=1:6 %back
    h=patch(x(:,i),y(:,i),z(:,i),'b');
    %set(h,'edgecolor','k');
    h=patch(xdx(:,i),ydy(:,i),zdz(:,i),'r');
   
    %set(h,'edgecolor','k');
end
end
alpha(.2);
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z');
Legend{1}='Undeformed';
Legend{2}='Deformed';
Legend{3}=strcat('max dx=',num2str(max(dx/scale)));
Legend{4}=strcat('max dy=',num2str(max(dy/scale)));
Legend{5}=strcat('max dz=',num2str(max(dz/scale)));
legend(Legend)

%legend('Original Nodes','Undeformed','Deformed')