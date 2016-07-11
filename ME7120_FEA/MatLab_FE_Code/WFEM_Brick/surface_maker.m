
global total_elements
global elnum
global nodes
global element
global surfaces

front=element(elnum).nodes(1:4);
back=element(elnum).nodes(5:8);
top=element(elnum).nodes([3 4 8 7]);
bottom=element(elnum).nodes([1 5 6 2]);
left=element(elnum).nodes([2 3 7 6]);
right=element(elnum).nodes([1 5 8 4]);

element(elnum).surfacenodes=[front;back;top;bottom;left;right];

total_elements=elnum;