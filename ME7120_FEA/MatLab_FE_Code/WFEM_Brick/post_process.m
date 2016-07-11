% post_process
global displacements
Xfull = full(X);
Xfull = Xfull(1:(length(nodes)*6));
displacements = reshape(Xfull,6,(length(nodes)));
displacements = displacements(1:3,:);

cubetry