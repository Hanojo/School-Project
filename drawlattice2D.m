function drawlattice2D(NODES,ELEMS)

% DRAWLATTICE2D(NODES,ELEMS) shows the lattice where the position of nodes
%         are given by the matrix NODES and the elements connecting the nodes
%         are given by the matrix ELEMS

Xnode=NODES(:,1);
Ynode=NODES(:,2);
figure
plot(Xnode,Ynode,'bo');
hold on;
for numelem=1:length(ELEMS)
   nodeleft=ELEMS(numelem,1);
   noderight=ELEMS(numelem,2);
   plot([Xnode(nodeleft) Xnode(noderight)],...
      [Ynode(nodeleft) Ynode(noderight)],'r');
end
hold off;
