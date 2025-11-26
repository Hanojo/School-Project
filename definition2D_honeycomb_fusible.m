function [NODES,ELEMS]=definition2D_honeycomb_fusible(L,W)

% fonction qui permet de definir les deux matrices NOEUD et ELEM qui vont
% caractériser un réseau initial en nid d'abeille de fusible de largeur L et de hauteur W 
%
% NODES est une matrice [NBOFNOEUD x 3] ou NBOFNOEUD est le nombre de noeuds
%       la première colonne NOEUD(:,1) donne la position X
%       la deuxième colonne NOEUD(:,2) donne la position Y
%       la troisième colonne NOEUD(:,3)donne la tension du noeud 
%                                       (nan à l'initiation)
%
% ELEM est une matrice [NBOFELEM x 5] ou NBOFELEM est le nombre de
%      resistance
%       la première colonne ELEM(:,1) donne le numéro du noeud en bas du
%                   résistance
%       la deuxième colonne ELEM(:,2) donne le numéro du noeud en haut du
%                   résistance
%       la troisième colonne ELEM(:,3) donne le courant passant 
%                                       (nan à l'initiation)
%       la quatrième colonne ELEM(:,4) donne le courant de claquage
%                                       (nan à l'initiation)
%       la cinquième colonne ELEM(:,5) donne la conductance
%                                       (1 à l'initiation)

% definition de NODES

Wd2=floor(L/2);
EPS=1e-4;

dY1=nan*ones(1,W);dY1(1:2:end)=0.5;dY1(2:2:end)=1;
dY1=[0,dY1(1:end-1)];

Y1=cumsum(dY1);

X1=zeros(size(Y1));
I=find(abs(round(Y1)-Y1)<EPS);
X1(I)=sqrt(3)/2;

Xnode=[];
Ynode=[];
for ind=0:Wd2-1
    Ynode=[Ynode,Y1];
    Xnode=[Xnode,X1+sqrt(3)*ind];
end
I=find(abs(round(Y1)-Y1)>EPS);
Ynode=[Ynode,Y1(I)];
Xnode=[Xnode,X1(I)+sqrt(3)*Wd2];

NODES=[Xnode(:) Ynode(:)];  

NODES=[NODES,nan*ones(size(NODES,1),1)];

% definition de ELEMS

ELEMS=[];
for ind=1:length(Xnode)
   Xcurrent=Xnode(ind);  % 
   Ycurrent=Ynode(ind);
   I=find( (abs((Xnode-Xcurrent).^2+(Ynode-Ycurrent).^2-1)<EPS) ...
       & (Ynode>=Ycurrent));
   ELEMS=[ELEMS;[ind*ones(length(I),1) I(:)]];
   %I=find( (abs((Xnode-Xcurrent).^2+(Ynode-Ycurrent).^2-1)<EPS) ...
   %   & (Xnode==Xcurrent)   (Ynode>=Ycurrent));
end

% On rajoute une 3eme colonne avec  le courant passant initialement à nan

ELEMS(:,3)=nan*ones(size(ELEMS,1),1);

% On rajoute une 4eme colonne avec le courant de claquage initialement à nan

ELEMS(:,4)=nan*ones(size(ELEMS,1),1);

% On rajoute une 5eme colonne avec la résistance initiale de 1

ELEMS(:,5)=ones(size(ELEMS,1),1);

% Origine au milieu du réseau

Xm=(min(NODES(:,1))+max(NODES(:,1)))/2;
Ym=(min(NODES(:,2))+max(NODES(:,2)))/2;

NODES(:,1)=NODES(:,1)-Xm;
NODES(:,2)=NODES(:,2)-Ym;
