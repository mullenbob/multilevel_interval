
function  [nn,nnA,ne,xA,yA,ndof,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1]=readtruss2(name,inp,fid)
%reading nodal information
[mat1] = fscanf(inp,'%d',2);
nnA = mat1(1); %number of nodes - assembled model
ne  = mat1(2); %number of elements
%x and y coodinates of nodes for assembled model
%xA = zeros(nnA,1);
%yA = zeros(nnA,1);
nn = 2*ne;    %number of nodes - EBE model
%x and y coodinates of nodes for EBE model
%xEBE = zeros(nn,1);
%yEBE = zeros(nn,1);
ndofA=2*nnA; %dof for assembled model
ndof =nn; %dof for EBE model  changed to just one dof per node 
%reading nodal information for the assembled model
[mat1] = fscanf(inp,'%d %f %f %d %d %e %e  %f %f',[9,nnA]);
mat1 = mat1';
xA = mat1(:,2);
yA = mat1(:,3);
%nodal information for assembled model
resxA = mat1(:,4); %x restraint for assembled model
resyA = mat1(:,5); %x restraint for assembled model
%x and y components of forces for assembled model
fxA = mat1(:,6);
fyA = mat1(:,7);
%uncertainty in forces for assembled model
betaxA=mat1(:,8);
betayA=mat1(:,9);

%A       = zeros(ne,1);
%E       = zeros(ne,1);
%Alfa    = zeros(ne,1);
%element = zeros(ne,2); %nodal connectivity matrix for EBE model
%--------------------------------------------------------------------------
%reading element information for the assembled model
[mat2] = fscanf(inp,'%d %d %d %f %e %f %d',[7,ne]);
mat2 = mat2';
elementA = mat2(:,2:3); %nodal connectivity matrix for assembled model
A       = mat2(:,4);
E       = mat2(:,5);
alfaA    = mat2(:,6); %uncertainty of E in assembled model
%--------------------------------------------------------------------------

gn      = mat2(:,7);  %group number
[mat3] = fscanf(inp,'%d',1);
ng = mat3(1);%Number of groups
if (ng >0)
[mat3] = fscanf(inp,'%d %f',[2,ng]);
mat3 = mat3';
gamma1 = mat3(:,2);% group uncertainty
% build new values of gammao using group data
gammao=zeros(ne,1);
for i=1:ne
    j=gn(i);
    if (j> 0)
        gammao(i)=gamma1(j);
    end
end
Gamma = infsup(1-gammao,1+gammao);
else
    gammao=zeros(ne,1);
    Gamma = infsup(1-gammao,1+gammao);
end
Alfa = infsup(1-alfaA,1+alfaA);
BetaxA = infsup(1-betaxA,1+betaxA);
BetayA = infsup(1-betayA,1+betayA);
% Print out input data now 
fprintf(fid,'File name %s \n',name);
fprintf(fid,'Nodal information - Assembled  model\n');
fprintf(fid,'Node   X      Y     Restraints       Fx          Fy           beta-x    beta-y\n');
for i=1:nnA
     fprintf(fid,'%2d   %4.1f   %4.1f    %d    %d     %9.1f    %9.1f           %5.3f      %5.3f\n',i,xA(i),yA(i),resxA(i),resyA(i),fxA(i),fyA(i),betaxA(i),betayA(i)); 
end
fprintf(fid,'Element information- Assembled model\n');
fprintf(fid,'Element    Nodes    Area            E        alfa      group no.    group gamma \n');
for i=1:ne
    fprintf(fid,'%2d        %2d  %2d    %f   %.3e   %7.5f         %2d           %7.5f\n',i,elementA(i,1),elementA(i,2),A(i),E(i),alfaA(i),gn(i), gammao(i));
end
fprintf(fid,'Groups Interval Solution\n');
fprintf(fid,'Number of Groups\n');
fprintf(fid,'%2d\n',ng);
fprintf(fid,'Group Number        Group Uncertainty   \n');
for i=1:ng
    fprintf(fid,'%2d                       %7.5f    \n',i,gamma1(i)); 
end
return
end