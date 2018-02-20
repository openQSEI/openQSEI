%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function: convert an image (unicorn) to a triangularized %%%
%%% discretization                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Danny Smyl
%%% Date: 10.2.2018
%%% Location: Aalto University, Espoo, Finland
%%% Corresponding Email: danny.smyl@aalto.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

img = imread('Unicorn.png');
A = img;
idx = find(A > 0);
A(idx) = 1;

[x,y]=meshgrid(1:size(A,1), 1:size(A,2));
x = 10*x/max(max(x));
y = 10*y/max(max(y));

x = reshape(x,[numel(x),1]);
y = reshape(y,[numel(y),1]);
A = reshape(A,[numel(A),1]);

F=scatteredInterpolant(x,y,double(A),'natural');
int=F(Ex,Ey);
Eunicorn = int;


save Eunicorn.mat Eunicorn

figure(1)
trisurf(Tri(IO,:),Tri.Points(:,1), Tri.Points(:,2),int),view(2),colormap('jet'),set(0,'defaulttextInterpreter','latex'),daspect([1 1 1]),shading interp, colorbar('eastoutside');caxis([cmin cmax]),axis tight,title('$E_{true}$','FontSize',14), axis off, box on;
