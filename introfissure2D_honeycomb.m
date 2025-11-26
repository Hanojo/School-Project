function [ELEMSnew,X0new,Y0new]=introfissure2D_honeycomb(NODES,ELEMS,X0,Y0)

% fonction qui introduit une fissure droite horizontale de pointe X0,Y0 
% en coupant les éléments adéquats
%
% NODES et ELEMS sont les matrices Noeuds et Eléments crées par les fonctions
% DEFINITION2D_##_###.m
%
% X0,Y0 est la position de la pointe de fissure
%
% ELEMSRES donne la matrice Eléments ainsi obtenue

EPS=1e-6;

ELEMSnew=ELEMS;
Xelem=0.5*(NODES(ELEMS(:,1),1)+NODES(ELEMS(:,2),1));
Yelem=0.5*(NODES(ELEMS(:,1),2)+NODES(ELEMS(:,2),2));


numelemref=(1:size(ELEMS,1))';
Xelemsnodes=NODES(ELEMS(:,2),1)-NODES(ELEMS(:,1),1);

% On identifie la pointe de fissure vraie,

I=find(abs(Xelemsnodes)<EPS);
Xelempossible=Xelem(I);
Yelempossible=Yelem(I);
[temp,I]=min((Xelempossible-X0).^2+(Yelempossible-Y0).^2);
X0new=Xelempossible(I);
Y0new=Yelempossible(I);

% On identifie les éléments à casser

indelem2withdraw=find((abs(Yelem-Y0new)<EPS) & (Xelem<X0new));
numelem2withdraw=numelemref(indelem2withdraw);

% On garde uniquement les autres dans ELEMSRES

numelemremaining=setdiff(numelemref,numelem2withdraw);
ELEMSnew=ELEMS(numelemremaining,:);