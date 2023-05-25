clear
clearvars -global

global tx;
global ty;
global force;
global fxA;
global fyA;
global iii;
global nnd
global nel
global ndof
global telement
global tE;
global tA;
global resxA;
global resyA;
global delta;
global iii;
global ngroups;
global gamma;
global groupid;
global alphaA;
global betax;
global betay;
global numintv;
global nintvload;
global nintvgroup;
global groupid;
global nodeid;
global beta;
format long;
%set counter for evaluation of function
iii=0;
name="popova5group"

inp= fopen(name+'.inp','r');
out =fopen(name+'optga.out','w');


readtrussx(inp,out);
lb=zeros(nel+nintvgroup+nintvload,1);
ub=lb;
x0=lb;

 %set initial values for modulus values
for i1=1:nel
    x0(i1)=tE(i1)/1.E+6;
    lb(i1)=x0(i1)*(1.-alphaA(i1)); 
    ub(i1)=x0(i1)*(1.+ alphaA(i1));
end
if (nintvgroup > 0)
    for i1=1:nintvgroup
    i2=i1+nel;
    x0(i2)=1.;
    lb(i2)=(1-gamma(i1));
    ub(i2)=(1.+gamma(i1));
    end
end
if (nintvload > 0)
for i1=1:nintvload
    i2=i1+nel+nintvgroup;
    x0(i2)=1.;
    lb(i2)=(1-beta(i1));
    ub(i2)=(1.+beta(i1));
end
end
tic
A=[];
Aeq=[];
b=[];
beq=[];
disp=trusssolve(x0)/10000.  % centered value check
tic
%options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxFunctionEvaluations',30000000);
%options = optimoptions('fmincon','PlotFcn','optimplotconstrviolation','Algorithm','sqp','Display','iter','MaxFunctionEvaluations',30000000);
%[x,fval,exitflag,output] = fmincon(@funx,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);
nvars=nel+nintvload+nintvgroup

options = optimoptions('ga','ConstraintTolerance',1e-6,'FunctionTolerance',1.e-10);
[x,fval,exitflag,output,population,scores] = ga(@funx,nvars,A,b,Aeq,beq,lb,ub,@nonlcon,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
x(nvars)
fprintf(out,' ga solution min %15.10e max %15.10e  fvalue min %15.10e  time %s function evals %d\n', min(x(1:nel)),max(x(1:nel)),fval/10000.,toc,iii);
iii=0;
tic;
[x,fval,exitflag,output,population,scores] = ga(@funy,nvars,A,b,Aeq,beq,lb,ub,@nonlcon,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' ga solution min %15.10e max %15.10e  fvalue max %15.10e  time %s function evals %d\n', min(x(1:nel)),max(x(1:nel)),-fval/10000.,toc,iii);
fprintf(out,'\nNodal information - optimization model\n');
fprintf(out,'Node   X      Y     Restraints       Fx          Fy           U-x    U-y\n');
ii=1;

for i=1:nnd
    
     fprintf(out,'%2d   %4.1f   %4.1f    %d    %d     %9.1f    %9.1f           %10.8g      %10.8g\n',i,tx(i),ty(i),resxA(i),resyA(i),fxA(i),fyA(i),delta(ii),delta(ii+1)); 
ii=ii+2;
end


    function [c,ceq] = nonlcon(x);
c = [];
ceq = [ ];
end
function disp=funx(optin);
global iii;
iii=iii+1;
% if (mod(iii,1000)== 0)
%    
% min(optin)
% max(optin)
% end
disp=trusssolve(optin);
return
end
function disp=funy(optin);
global iii;
iii=iii+1;

disp=-trusssolve(optin);
return
end
function readtrussx(inp,out)
   
% set up global variables to exchange inifor between functions

global tx;
global ty;
global force;
global nnd
global nel
global ndof
global telement
global tE;
global tA;
global fxA;
global fyA;
global resxA;
global resyA;
global ngroups;
global gamma;
global groupid;
global alphaA;
global betax;
global betay;
global numintv;
global nintvload;
global nintvgroup;
global groupid;
global nodeid;
global beta;



% read in data in same forame as 
                                                                            
[mat1] = fscanf(inp,'%d',2);
nnd= mat1(1) %number of nodes - assembled model
nel  = mat1(2) %number of elements
ndof=2*nnd

%x and y coodinates of nodes for assembled model
tx = zeros(nnd,1);
ty = zeros(nnd,1);
resxA=zeros(nnd,1);
resyA=zeros(nnd,1);
fxA=zeros(nnd,1);
fyA=zeros(nnd,1);

%reading nodal information for the assembled model
[mat1] = fscanf(inp,'%d %f %f %d %d %e %e  %f %f',[9,nnd]);
mat1 = mat1';
tx = mat1(:,2);
ty = mat1(:,3);
%nodal information for assembled model
resxA = mat1(:,4); %x restraint for assembled model
resyA = mat1(:,5); %x restraint for assembled model
%x and y components of forces for assembled model
fxA = mat1(:,6);
fyA = mat1(:,7);
betax = mat1(:,8);
betay = mat1(:,9);
% count number of interval loads
nintvload=0;
for i1=1:nnd
    if (betax(i1) > 0)
        nintvload=nintvload+1;
        nodeid(nintvload)=i1*2-1;
        beta(nintvload)=betax(i1);
    end
    if (betay(i1) > 0) 
        nintvload=nintvload+1;
        nintvload=nintvload+1;
        nodeid(nintvload)=i1*2;
        beta(nintvload)=betay(i1);
    end
end
%reading element information for the assembled model
tA       = zeros(nel,1);
tE       = zeros(nel,1);
telement = zeros(nel,2);
[mat2] = fscanf(inp,'%d %d %d %f %e %f %d',[7,nel]);
mat2 = mat2';
telement = mat2(:,2:3); %nodal connectivity matrix for assembled model
tA       = mat2(:,4);
tE       = mat2(:,5);
alphaA    = mat2(:,6); %uncertainty of E in assembled mode
groupid = mat2(:,7); 
% read in group information
ngroups=fscanf(inp,'%d',1);

nintvgroup=0;
if (ngroups > 0)
  [mat4,kode] = fscanf(inp,'%d %f',[2,ngroups]);

mat4=mat4';
gamma= mat4(:,2);

    for i1=1:ngroups
    if (gamma(i1) > 0)
        nintvgroup=nintvgroup+1;
    end  
    end
end
return
end
function temp=trusssolve(optvar)
global ndof;
global delta;
global nel;
global tx;
global ty;
global force;
global nnd
global nel
global ndof
global telement
global tE;
global tA;
global fxA;
global fyA;
global resxA;
global resyA;
global ngroups;
global gamma;
global groupid;
global alphaA;
global betax;
global betay;
global numintv;
global nintvload;
global nintvgroup;
global groupid;
global nodeid;
global beta;

K=zeros(ndof);
  %assenble stiffness matrix
  for e=1:nel
      tEnew=optvar(e);
      ig=groupid(e);
      if (nintvgroup > 0 && ig > 0 && optvar(nel+ig)>0)
          tEnew=tEnew*optvar(nel+ig);
      end
     i=telement(e,1);
     j=telement(e,2);
     txx=tx(j)-tx(i);
     tyy=ty(j)-ty(i);
     l=sqrt(txx^2+tyy^2);
     c=txx/l;
     s=tyy/l;
     k11=(tEnew*1.0e6*tA(e)/l)*[(c)^2 s*c;s*c s^2];
    K((2*i-1):2*i,(2*i-1):2*i)=K((2*i-1):2*i,(2*i-1):2*i)+k11;
    K((2*j-1):2*j,(2*j-1):2*j)=K((2*j-1):2*j,(2*j-1):2*j)+k11;
    K((2*i-1):2*i,(2*j-1):2*j)=K((2*i-1):2*i,(2*j-1):2*j)-k11;
    K((2*j-1):2*j,(2*i-1):2*i)=K((2*j-1):2*j,(2*i-1):2*i)-k11;
  end
%
% force fixed dof
for e=1:nnd
    
if (resxA(e) == 0)
        ii=e*2-1;
        for i=1:nnd
            i1=i*2-1;
            K(ii,i1)=0.;
            K(i1,ii)=0.;
            i1=i1+1;
            K(ii,i1)=0.;
            K(i1,ii)=0.;
        end
        K(ii,ii)=1.0;
end
if (resyA(e) == 0)
        ii=e*2;
        for i=1:nnd
            i1=i*2-1;
            K(ii,i1)=0.;
            K(i1,ii)=0.;
            i1=i1+1;
            K(ii,i1)=0.;
            K(i1,ii)=0.;
            end
        K(ii,ii)=1.0;
end
end
% build force vector
ii=0;
for i=1:nnd
ii=ii+1;
force(ii)=fxA(i);
ii=ii+1;
force(ii)=fyA(i);
end
if (nintvload > 0)
for i=1:nintvload
     i2=i+nel+nintvgroup;
   
    ii=nodeid(i);
    force(ii)=force(ii)*optvar(i2);
end
end
delta=K\force';
temp=delta(2*nnd-1)*10000;
return  
end
