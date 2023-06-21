%LINEAR INTERVAL FINTE ELEMaANALYSIS OF A Three dimensional  TRUSS
%SECONDARY VARIABLES WITH THE SAME ACCURACY OF THE PRIMARY ONES
%A PROGRAM BY Dr. RAMA RAO MALLELA , Prof. RAFI MUHANNA AND Prof.  ROBERT MULLEN
% addpath 'G:\My Drive\Documents\MATLAB\Intlab_V8\Intlab_V8'
%  startintlab();
 intvalinit('DisplayInfsup');
clear
clc
format long
tic
intvalinit('DisplayInfsup');
name="threeD1"
inp= fopen(name+'.inp','r');
out =fopen(name+'.out','w');

fid=out;
%reading nodal information
[mat1] = fscanf(inp,'%d',2);
nnA = mat1(1); %number of nodes - assembled model
ne  = mat1(2); %number of elements
%x and y coodinates of nodes for assembled model
xA = zeros(nnA,1);
yA = zeros(nnA,1);
zA = zeros(nnA,1);
nn = 2*ne;    %number of nodes - EBE model
%x and y coodinates of nodes for EBE model
xEBE = zeros(nn,1);
yEBE = zeros(nn,1);
zEBE = zeros(nn,1);
ndofA=3*nnA; %dof for assembled model
ndof =3*nn; %dof for EBE model
%reading nodal information for the assembled model
[mat1] = fscanf(inp,'%d %f %f %f %d %d %d %e %e %e  %f %f %f',[13,nnA]);
mat1 = mat1';
xA = mat1(:,2);
yA = mat1(:,3);
zA = mat1(:,4);
%nodal information for assembled model
resxA = mat1(:,5); %x restraint for assembled model
resyA = mat1(:,6); %y restraint for assembled model
reszA = mat1(:,7); %z restraint for assembled model
%components of forces for assembled model
fxA = mat1(:,8);
fyA = mat1(:,9);
fzA = mat1(:,10);
%uncertainty in forces for assembled model
betaxA=mat1(:,11);
betayA=mat1(:,12);
betazA=mat1(:,13);

bnA = 0; %number of restrained dof in assembled model
for i=1:nnA
   if(resxA(i)==0) %x-dof is restrained
        bnA = bnA+1;
        ifixA(bnA) = 3*i-2;
    end
   if(resyA(i)==0) %y-dof is restrained
        bnA = bnA+1;
        ifixA(bnA) = 3*i-1;
   end
   if(reszA(i)==0) %z-dof is restrained
        bnA = bnA+1;
        ifixA(bnA) = 3*i;
   end
end
A       = zeros(ne,1);
E       = zeros(ne,1);
Alfa    = zeros(ne,1);
element = zeros(ne,2); %nodal connectivity matrix for EBE model
%--------------------------------------------------------------------------
%reading element information for the assembled model
[mat2] = fscanf(inp,'%d %d %d %f %e %f %d',[7,ne]);
mat2 = mat2';
elementA = mat2(:,2:3); %nodal connectivity matrix for assembled model
A       = mat2(:,4);
E       = mat2(:,5);
alfaA    = mat2(:,6); %uncertainty of E in assembled model
%--------------------------------------------------------------------------
Del_alpha = infsup(-alfaA,alfaA);
gn      = mat2(:,7);  %group number
% gammao  = mat2(:,8);   % group uncertainty, i will work to remove these
% in the future.  done!  RLM 9/6/22
[mat3] = fscanf(inp,'%d',1);
ng = mat3(1);%Number of groups
if (ng >0)
[mat3] = fscanf(inp,'%d %f',[2,ng]);
mat3 = mat3';
gamma1 = mat3(:,2);% group uncertainty
% build new values of gammao using group data  
% conversion
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
%start timing
tic
%--------------------------------------------------------------------------
Alfa = infsup(1-alfaA,1+alfaA);
BetaxA = infsup(1-betaxA,1+betaxA);
BetayA = infsup(1-betayA,1+betayA);
BetazA = infsup(1-betazA,1+betazA);
nforce=0;
ForceA=infsup(zeros(3*nnA,1),zeros(3*nnA,1)); %allocate space for ForceA  RLM 9/8/22
for i=1:nnA
  if(abs(fxA(i))>0.0)
      nforce = nforce+1;
  end
  if(abs(fyA(i))>0.0)
      nforce = nforce+1;
  end
if(abs(fzA(i))>0.0)
      nforce = nforce+1;
  end
  ForceA(3*i-2) = BetaxA(i)*fxA(i);
  ForceA(3*i-1)   = BetayA(i)*fyA(i);
  ForceA(3*i)   = BetazA(i)*fzA(i);
end
nE = 0; %number of elements with uncertain E
for i=1:ne
    if(alfaA(i)>0.0)
        nE = nE+1;
    end
end
for e=1:ne
    element(e,1)=2*e-1;
    element(e,2)=2*e;
end
fx = zeros(nn,1);
fy = zeros(nn,1);
fz = zeros(nn,1);
loadflag = zeros(nnA,1);
for i=1:nnA
    loadflag(i) = 0;
    if((fxA(i)~=0)||(fyA(i)~=0))||fzA(i) ~=0
        loadflag(i) =1;
    end
end
ebenode=zeros(2*ne,1);
resxEBE=zeros(2*ne,1);
resyEBE=zeros(2*ne,1);
reszEBE=zeros(2*ne,1);
for e=1:ne
   i1 = elementA(e,1); %nodes for assembled model
   i2 = elementA(e,2); 
   j1 = 2*e-1;%nodes for EBE model
   j2 = 2*e;
   ebenode(j1)=i1;
   ebenode(j2)=i2;
   xEBE(2*e-1) = xA(i1);
   yEBE(2*e-1) = yA(i1);
   zEBE(2*e-1) = zA(i1);
   xEBE(2*e)   = xA(i2);
   yEBE(2*e)   = yA(i2);
   zEBE(2*e) =   zA(i2);
   %restraints for EBE model
   resxEBE(2*e-1) = resxA(i1);
   resyEBE(2*e-1) = resyA(i1);
   reszEBE(2*e-1) = reszA(i1);
   resxEBE(2*e)   = resxA(i2);
   resyEBE(2*e)   = resyA(i2);
   reszEBE(2*e)   = reszA(i1);
   %forces for EBE model
   if(loadflag(i1)==1)
       fx(j1)=fxA(i1);
       fy(j1)=fyA(i1);
       fz(j1)=fzA(i1);
       loadflag(i1)=0;
   end
   if(loadflag(i2)==1)
      fx(j2)=fxA(i2);
      fy(j2)=fyA(i2); 
      fz(j2)=fzA(i2);
      loadflag(i2)=0;
   end
end
betaxEBE = zeros(nn,1);
betayEBE = zeros(nn,1);
betazEBE = zeros(nn,1);
%load uncertainty for EBE model
for i=1:nn
  j=ebenode(i);
  betaxEBE(i) = betaxA(j);
  betayEBE(i) = betayA(j);
  betazEBE(i) = betazA(j);
end
%determining the list of coincident nodes for new approach
ic = 0;
for i=1:nn
    code = 0;
    xi = xEBE(i); %coordinates of i'th node
    yi = yEBE(i); 
    zi = zEBE(i);
    for j=1:nn
              xj = xEBE(j); %coordinates of j'th node
              yj = yEBE(j); 
               zj = zEBE(j);
              if (abs(xi-xj)<1.0e-6)&&(abs(yi-yj)<1.0e-6)&&(abs(zi-zj)<1.0e-6)
                %for coincident nodes which are not hinged nodes
                if(code==0)
                    ic = ic+1;
                    icvec(ic,1)=j;
                    icvec(ic,2)=i;
                    code =1; %do only first one found
                end
              end
    end
end
icvec
k=0;
codevec = zeros(ic,1);
for i=1:ic
    if (codevec(i)==0)
        k=k+1;
        i1= icvec(i,1);
        i2= icvec(i,2);
        if(i1~=i2)  
          list(k,1) = i1;
          list(k,2) = i2;
          jj=2;
        else
          list(k,1) = i1;
          jj=1;
        end
        for j =i+1:ic %search the remaining rows
            if(i1==icvec(j,1)) %common node
                 jj=jj+1;
                 list(k,jj) = icvec(j,2);
                 i2= icvec(j,2);
                  codevec(j)=1;
             end
        end            
    end
end
list
codevec
[nf i2]  = size(list)
codefree = ones(nf,3); % BC at common nodes,unrestrained by default
for i=1:nf  %loop on the total number of free nodes
   j1 = ebenode(list(i,1));
   if(resxA(j1)==0)
        codefree(i,1)=0; %free node is restrained along X direction
   end
   if(resyA(j1)==0)
        codefree(i,2)=0; %free node is restrained along Y direction
   end
   if(reszA(j1)==0)
        codefree(i,3)=0; %free node is restrained along Z direction
   end
end
codefree;
for i=1:nf
    j1=0;
    for j=1:i2
        lnode = list(i,j);
        if(lnode==0)
            elemlist(i,j)=0;
        end
        for e=1:ne
           if((lnode==(2*e-1))||(lnode==(2*e))) 
                   j1 = j1+1;
                   elemlist(i,j1) = e;
           end
        end
    end
end
elemlist;
ndof1 = 3*nf; %global dof for free nodes  % check 3 or 2
nzeros = 0;
for i=1:nf
    for j=1:i2
        if(list(i,j)==0)
            nzeros= nzeros+1;
        end
    end
end
nlamda = nf*i2-nzeros;
ndofnew = ndof+ndof1+nlamda+2*ne;
fprintf(1,'ndof = %d ndof1=%d nlamda= %d ndofnew=%d\n',ndof,ndof1,nlamda,ndofnew);
a = zeros(ne,ndofnew);
aA = zeros(ne,ndofA);
for e=1:ne
  conn=element(e,:);
  connA=elementA(e,:);
  x1 = xEBE(conn(1));
  y1 = yEBE(conn(1));
  z1 = zEBE(conn(1));
  x2 = xEBE(conn(2));
  y2 = yEBE(conn(2));
  z2 = zEBE(conn(2));
  le=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
  Le(e) = le;
  Cx=(x2-x1)/le;
  Cy=(y2-y1)/le;
  Cz=(z2-z1)/le;
  cosine(e) = Cx;
  sine(e) = Cy;
  sine2(e) = Cz;    %see why needed
  
  dof=[3*conn(1)-2 3*conn(1)-1 3*conn(1) 3*conn(2)-2 3*conn(2)-1 3*conn(2)];
  dofA=[3*connA(1)-2 3*connA(1)-1 3*connA(1) 3*connA(2)-2 3*connA(2)-1 3*connA(2)];
   %____________global transformation matrix_________________
  % R=[ c -s  0  0;
  %     s  c  0  0;
  %     0  0  c -s;
  %     0  0  s  c];
  % a(e,dof)= eye(4)*[1 0 -1  0]'; % RLM I do not see why we need the identy matrix
  % aA(e,dofA)= R*[1 0 -1  0]'; % again   why not just vector [c s -c -c]'
   a(e,dof) = [ 1 0 0 -1 0 0];
  aA(e,dofA)= [Cx Cy Cz -Cx -Cy -Cz];
  %____________diagonal member matrix_______________________
  % elemiminate ki and only use Ki  RLM 5/26/2023
  %Ki(e,e)=diag(A(e)*(E(e)*Alfa(e))/le);
 % Ko = mid(Ki);
  %obtaining global stiffness matrix using tranformation matrix  RLM should
  %remove the ,0 inthe diag calls as it defaults to being on the diagonal
 Doo(e,e)=diag([A(e)*E(e)/le]); 
 Do(e,e)=diag([A(e)*E(e)/le]*infsup(1.0,1.0));
 D(e,e)=diag([(A(e)*E(e)/le)*Alfa(e)]);   %I do not see any difference between D and Ki ?? RLM 5/26/2023
  % I replaced all Ki with D RLM 5/26/23
 DoD(e,e)=Do(e,e)-D(e,e) ; 
 Dom(e,e) = diag([(E(e)*A(e)/le)*(1+(alfaA(e)*gammao(e)))],0);  
end

%----------------formulation of the new approach-------------------------
[nf i2]  = size(list);
for i=1:nf
    list1 = list(i,:);
    list1 = list1';
    for j=1:size(list1)
      j1= list1(j);
      forceEBE(3*i-2)  = infsup(1-betaxEBE(j1),1+betaxEBE(j1))*fx(j1);
      forceEBE(3*i-1)    = infsup(1-betayEBE(j1),1+betayEBE(j1))*fy(j1);
      forceEBE(3*i  )    = infsup(1-betazEBE(j1),1+betazEBE(j1))*fz(j1);
      break;
    end
end
for i=1:ndofnew
    Force(i) = infsup(0.0,0.0);
end
for i=ndof+1:ndof+ndof1
 Force(i) = forceEBE(i-ndof);
end
for i=ndof+ndof1+1:ndofnew
    Force(i)=infsup(0.0,0.0);
end
Force=Force';
cons = zeros(2*ne+nlamda,ndofnew);%*infsup(0,0);
%filling in the entries of constraint matrix
%this matrix relates the displacements of coincident nodes along local axes
%with the global displacements of each of the corresponding free nodes
ilamda = 0;
[nrows ncols] = size(list);
for i=1:nrows
    for j=1:ncols
        lnode = list(i,j); %pick up a node from the list of coincident nodes
        for e=1:ne
               %determine to which element does this node belong to
               if(lnode==(2*e-1))||(lnode==(2*e)) 
                  %determining the far node of the element for which the
                  %near node is a coincident node 
                  if (lnode==(2*e-1))
                      farnode = 2*e;
                  end
                  if (lnode==(2*e))
                      farnode = 2*e-1;
                  end
                  ilamda = ilamda+1;
                  cons(ilamda,3*lnode-2) = 1.0;
                  cons(ilamda,ndof+3*i-2)= -cosine(e); 
                  cons(ilamda,ndof+3*i-1)  = -sine(e); 
                  cons(ilamda,ndof+3*i)  = -sine2(e); 
               end
        end  
    end
end
for e=1:ne
    node1 = 2*e-1;
    node2 = 2*e;  
    cons(nlamda+e,2*node1-1) = -1/Le(e);
    cons(nlamda+e,2*node2-1) =  1/Le(e);
end
%including Identity matrix of size ne
cons(nlamda+ne+1:nlamda+2*ne,ndofnew-2*ne+1:ndofnew-ne) = -eye(ne);
cmat = zeros(ndofnew,ndofnew);%*infsup(0,0);
for i=1:nlamda+2*ne
    %i
   for j=1:ndofnew
       %j
      cmat(ndof+ndof1+i,j)=cons(i,j);
      %cmat
   end
end
% make cmat symmetric ???
mat = cmat';
cmat = cmat+mat;
%-------------------------------------------------------------------
%application of boundary conditions
%-------------------------------------------------------------------
%NOTE : for the 4x4 local stiffness matrix of truss element, the dof
%corresponding to local y and axis are restrained
bn1 = 0;
resy2 = resyEBE;
resz2 = reszEBE;
for i=1:nn
    resyEBE(i) = 0;
    reszEBE(i) = 0;
end
for i=1:nn
   if(resyEBE(i)==0)
        bn1 = bn1+1;
        ifix2(bn1) = 3*i-1;
   end
   if(reszEBE(i)==0)
        bn1 = bn1+1;
        ifix2(bn1) = 3*i;
   end
end
resyEBE=resy2;
reszEBE=resz2;
%-------------------------------------------------------------------
%NOTE : now boundary conditions are imposed at the free nodes along global dof
bn2 = 0;
for i=1:nf
    j = ebenode(list(i,1));
    if(codefree(i,1)==0)  %restrained free node
         if(resxA(j)==0)
           bn2 = bn2+1;
           ifix3(bn2) = ndof+ 3*i-2;
         end
    end
    if(codefree(i,2)==0)  %restrained free node
        if(resyA(j)==0)
          bn2 = bn2+1;
          ifix3(bn2) = ndof+ 3*i-1;
        end
    end
    if(codefree(i,3)==0)  %restrained free node
        if(reszA(j)==0)
          bn2 = bn2+1;
          ifix3(bn2) = ndof+ 3*i;
        end
    end
end
%------------------------------------------------------------------------
%K = a'*Ki*a;  %replace Ki with D
K = a'*D*a;
K = K+cmat;
Km = a'*Dom*a;
Km = Km+mid(cmat);  %Changed to mid for Km to be non interval
% imposing restraint along local y and zaxes
K(ifix2,:)=zeros(bn1,ndofnew);
K(:,ifix2)=zeros(ndofnew,bn1);
K(ifix2,ifix2)=eye(length(ifix2));
% imposing BC along y axes for groups
Km(ifix2,:)=zeros(bn1,ndofnew);
Km(:,ifix2)=zeros(ndofnew,bn1);
Km(ifix2,ifix2)=eye(length(ifix2));
%imposing restraints on free nodes
K(ifix3,:)=zeros(bn2,ndofnew);
K(:,ifix3)=zeros(ndofnew,bn2);
K(ifix3,ifix3)=eye(length(ifix3)); 
% calculating the inverse of mid-point matrices K
C = inv(mid(K));
% imposing BC on free nodes
Km(ifix3,:)=zeros(bn2,ndofnew);
Km(:,ifix3)=zeros(ndofnew,bn2);
Km(ifix3,ifix3)=eye(length(ifix3));
% calculating the inverse of mid-point matrices Km
Cm =inv(Km);
%*********************   Initial Enclosure  *******************
w=ones(ne,1);
w1=w-mag(DoD)*mag(a*C*a')*w;
w2=mag(DoD)*mag(a*C*Force);
% allocate alpha for efficency
alpha=zeros(ne,1);
for i=1:ne
    alpha(i)=(w2(i)/w1(i));
end
alphamax=max(alpha(:,:));           
dd=infsup(-alphamax*w,alphamax*w);  
U1 = C*Force + C*a'*dd;
%************************************************************
for i=1:bn1    
  a(:,ifix2(i))=0; 
end 
v=a*U1;
d=(DoD)*v;
v(:,1)=v;
d(:,1)=d; 
for i=1:10
    v(:,i+1)=intersect(((a*C)*Force+(a*C*a')*d(:,i)),(v(:,i)));
    d(:,i+1)=intersect((DoD*v(:,i+1)),(d(:,i)));
end
u1=(C*Force)+(C*a')*d(:,11); %Du1 calculated using first approach-Element uncertainty only
%--------------------------------------------------------------------------
%---------------------------------------------------------------
%for large uncertainties RLM commented out  9/9/22 
% c=Doo*a*C*Force;
% for i=1:ne
%     sum=0;
%     for k=1:ne
%        sum=sum+mag(c(k))^2/inf(D(k,k));
%     end
% sum=sum*sup(D(i,i))+mag(c(i))^2*(1-sup(D(i,i))/inf(D(i,i)));
% sum=sqrt(sum);
% z(i)=1/2*(c(i)+infsup(-1,1)*sum);
% vnew(i)=z(i)/D(i,i);
% Dnew(i)=(Doo(i,i)/D(i,i)-1)*z(i);
% end
% v(:,1)=vnew;
% d(:,1)=Dnew; 
% %The intersection iteration
% for i=1:bn1    
%   a(:,ifix2(i))=0; 
% end 
% for i=1:10
%     v(:,i+1)=intersect(((a*C)*Force+(a*C*a')*d(:,i)),(v(:,i)));
%     d(:,i+1)=intersect((DoD*v(:,i+1)),(d(:,i)));
% end
%u1large=(C*Force)+(C*a')*d(:,11); %Du1 calculated when large uncertainty is present
% not sure why there are two solutions for u1  (large uncertainty?????)
% difference is 10^-22   I will be commentting this section out RLM 9/9/22
%--------------------------------------------------------------------------
% Group Uncertainty
if (ng > 0)
LambdaG=zeros(ne,ng);

    for i = 1:ne
    j = gn(i);
    if j~=0
        LambdaG(i,j)=1;
    end
end
%LambdaG
gamma = gamma1*(infsup(-1,1));
end
% nc =size(LambdaG,2)
% gamma=zeros(nc,1);
if ng == 0
    gamma=zeros(nc,1);
    Del_gamma = gamma;
else
for i= 1:ng
    Del_gamma(i)= gamma(i);
end
Del_gamma = Del_gamma';
Del_gammas=1-Del_gamma;
end
%Del_gamma = Del_gamma';
Del_gammas=1-Del_gamma;
LambdaE = eye(ne,ne);  %may not need if E is EBE but keep for reading equations
%--------------------------------------------------------------------------
vs = (a*C)*Force;
%vss = (a*C)*Force;
vg = (a*C)*Force;
% changed to Alfax to save old variable Alfa for combination solution RLM
% 9/8/22
Alfax = diag(DoD);
for i = 1:10
    %vss=hull(((a*Cm)*Force-((a*Cm*a'*Ko)*diag(vss)*LambdaE*Del_alpha)),vss);
    vg = hull(((a*C)*Force - (((a*C*a')*diag(vg))*LambdaE)*Alfax),vg);
   % replace Ko with Do  vs = hull(((a*Cm)*Force-((a*Cm*a'*Ko)*diag(vs)*LambdaE*Del_alpha)+((a*Cm*a'*Ko)*diag(vs)*LambdaG*Del_gamma)),vs);
    vs = hull(((a*Cm)*Force-((a*Cm*a'*Do)*diag(vs)*LambdaE*Del_alpha)+((a*Cm*a'*Do)*diag(vs)*LambdaG*Del_gamma)),vs);
end
% replace Ko with Do RLM 5/26/23
% u1s = (Cm*Force-((Cm*a'*Ko)*diag(vs)*LambdaE*Del_alpha)+((Cm*a'*Ko)*diag(vs)*LambdaG*Del_gamma));
u1s = (Cm*Force-((Cm*a'*Do)*diag(vs)*LambdaE*Del_alpha)+((Cm*a'*Do)*diag(vs)*LambdaG*Del_gamma));
%u1ss = (Cm*Force-((Cm*a'*Ko)*diag(vss)*LambdaE*Del_alpha));
ug = C*Force + C*a'*diag(vg)*LambdaE*Alfax; 
IFEAtime = toc
%--------------------------------------------------------------------------
%ug is the vector of global displacements at free nodes
for i=ndof+1:ndof+ndof1
    ug(i-ndof) = u1(i);
    usg(i-ndof) = u1s(i);
end
%constructing the overall displacement vector {Du}
for i=1:ndof
    u(i) = infsup(0.0,0.0);
    us(i)= infsup(0.0,0.0);
end
for i=1:nn
      %if the node is among the list of restrained nodes
      if(resxEBE(i)==0)||(resyEBE(i)==0)
          if(resxEBE(i)==1)
            u(2*i-1) = infsup(0.0,0.0);
            us(2*i-1) = infsup(0.0,0.0);
            
          end
          if(resyEBE(i)==0)
              u(2*i) = infsup(0.0,0.0);
              us(2*i) = infsup(0.0,0.0);
          end  
      end
      %if the node is among the list of coincident nodes 
          for i1=1:nf
              for i3=1:size(list,2)
                  if(i==list(i1,i3))
                      u(2*i-1) = ug(2*i1-1);
                      u(2*i) = ug(2*i1);
                      us(2*i-1) = usg(2*i1-1);
                      us(2*i) = usg(2*i1);
                      
                  end
              end
          end
end
for i=ndof+ndof1+1:ndof+ndof1+nlamda
    j=i-ndof-ndof1;
    lamda(j) =u1(i);   % this is the EBE only solution  
    lamdas(j) =u1s(i);  % this is the EBE + Groups solution
%     fprintf(1,'[%6.0f,%6.0f]\n',inf(lamda(j)),sup(lamda(j)));
end
%assigning contents of lamda vector to the element force vector
[nrows ncols] = size(list);
for i=1:e
    for j=1:2
        eforce(i,j) = infsup(0.0,0.0);
        eforces(i,j) = infsup(0.0,0.0);
    end
end
index =0;
for i=1:nrows
   for j=1:ncols
     lnode= list(i,j);
     if(lnode~=0)
         index =index+1;
          elem = elemlist(i,j);
          if(2*elem-1==lnode)
             eforce(elem,1)=lamda(index);
             eforces(elem,1)=lamdas(index);
          end
          if(2*elem==lnode)
             eforce(elem,2)=lamda(index);
             eforces(elem,2)=lamdas(index);
          end
     end
   end
end
% for i=1:ne
%     fprintf(1,'element=%d [%12.5f,%12.5f]  [%12.5f,%12.5f]\n',i,inf(eforce(i,1)),sup(eforce(i,1)),inf(eforce(i,2)),sup(eforce(i,2)));
% end
% %mat = [  1   0   -1    0
%          0   0    0    0
%         -1   0    1    0
%          0   0    0    0]; 
%locations of non-zero loads in assembled force vector
 j=0;
 intforce = 0;
 for i=1:ndof1
   
    if(mid(ForceA(i)~=0))
        j=j+1;
        locx(j)=i;
    end
    if((sup(ForceA(i))-inf(ForceA(i)))>0)
        intforce = intforce+1;
    end
 end
 ng1 = 0;
 for i = 1:ng
     if rad(gamma(i))>0
         ng1=ng1+1;
     end
 end
 fprintf(out,'fixed point interval solution time  %s\n',toc);
  %fprintf(fid,'Total combinations is %d \n',(2^ne)*(2^ng1)*(2^intforce));
 fprintf(1,'Total combinations is %d \n',(2^ne)*(2^ng1)*(2^intforce));
 fprintf(1,'Do you want to perform combinatorial solution?\n');
 option = input('Press 1 for Yes or 2 for No :');
  if (option==1) 
     tic;
jmax=(2^ne)*(2^ng1)*(2^intforce);

printfeq=1000;
if (jmax < 1000)
    printfeq=100;
end
if (jmax < 100)
    printfeq=10;
end
 fprintf(fid,'Combinatorial solution ran\n');
        %Combinatorial solution
     Q1=zeros(ne,2);
     Q2=zeros(ne,2);
     for i=1 :ne
         Q1(i,:)=[inf(E(i)*Alfa(i)),sup(E(i)*Alfa(i))];
     end
     %added multiplier for group term 9/6/22  RLM
     for i=1: ne
        Q2(i,:)=[inf(Gamma(i)),sup(Gamma(i))];
     end
     Q = zeros(ndof1,2);
     for i=1:ndof1
       Q(i,:)= [inf(ForceA(i)); sup(ForceA(i))];
     end
     %force combinations
     %using binary numbers to generate combinations of lower and upper bounds
     %of forces for loads =nforce  
     %rlm   this looks like a complex way to calculated 2^nforce-1
     bin1=ones(intforce,1);
     buf = int2str(bin1)';
     ncombF= bintodec(buf);
     bin2 = ones(nE,1);
     buf = int2str(bin2)';
     ncombE= bintodec(buf);
     %added bin3 for group data
     bin3 = ones(ng1,1);
     buf = int2str(bin3)';
     ncombG= bintodec(buf);
     fprintf(fid,'Number of combinations force %d, combinations E %d, combinations A (groups) %d\n',ncombF+1,ncombE+1,ncombG+1);
     if (ng1==0)
     ncombG=0;
     end

     if(nE==0)
         ncombE=0;
     end
     if(intforce==0)
         ncombF=0;
     end
     %generating combinations of upper and lower bounds of loads and E and
     %and groups 
     count=0;
%     countE=1;   %Not used   commented it out RLM 9/9/22
for ncE=0:ncombE
         combE =dectobin(ncE,ne);
         for i=1:ne
            if(combE(i)==0)
                E1(i)=Q1(i,1);
            else
                E1(i)=Q1(i,2);
            end
         end

for ncG=0:ncombG
             combG =dectobin(ncG,ng)';
             %modify area for group uncertainty RLM 9/8/22
            
         for i=1:ne
             ii=gn(i);
             A1(i)=A(i);
             if (ii > 0)
            if(combG(ii)==0)
                A1(i)=A(i)*Q2(i,1);
            else
                A1(i)=A(i)*Q2(i,2);
            end
             end
         end


for ncF=0:ncombF 
            count = count+1;
            combF =dectobin(ncF,nforce)';
            for i=1:ndof1
                P1(i)=0;
            end
            for i=1:nforce
                j=locx(i); %location of dof carrying non-zero load
               if(combF(i)==0)  %set load to lower bound of ForceA
                    P1(j)= inf(ForceA(j));
               else
                    P1(j)=sup(ForceA(j));%set load to upper bound of ForceA
               end
            end
            for e=1:ne
             conn1=elementA(e,:);
             dof2(:,:,e)=[3*conn1(1)-2 3*conn1(1)-1 3*conn1(1) 3*conn1(2)-2 3*conn1(2)-1 3*conn1(2)];
             x1=xA(conn1(1));
             y1=yA(conn1(1));
             z1=zA(conn1(1));
             x2=xA(conn1(2));
             y2=yA(conn1(2));
             z2=zA(conn1(2));
             le=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
             Cx=(x2-x1)/le;
             Cy=(y2-y1)/le; 
             Cz=(z2-z1)/le; 
             %changed T to a b matrix 
             T(:,:,e)= [-Cx  -Cy -Cz  Cx Cy Cz ];
                      
             % change to A1 for groups
             kc(e)=A1(e)*(E1(e))/le; 
            end
            Kic(:,:)=diag(kc(:));
            Kicn(:,:)=aA'*Kic(:,:)*aA;
            Kicn(ifixA,:)=0;                   
            Kicn(:,ifixA)=0;                    
            Kicn(ifixA,ifixA)=eye(length(ifixA));
%             ucvec =inv(Kicn)*P1';
             ucvec =Kicn\P1';
             % the stress and strain calculations are strange as it is
             % using a 4x4 matrix for one value (the third which is the
             % interanl force at the third dof in the local system. RLM
             % 9/9/22
            for e=1:ne
              ke1 = kc(e)*mat;  
              ucx = ucvec(dof2(:,:,e));  %changed to ucx to not repeat uc
              g1(e) = kc(e)*T(:,:,e)*ucx;
              %changed to use A1 to account for groups
              strain1(e) = g1(e)/(A1(e)*E1(e));
              if(count==1)
                  gmin(e)= g1(e);
                  gmax(e)= g1(e);
                  strainmin(e) = strain1(e);
                  strainmax(e) = strain1(e);
              else
                gmin(e) = min(g1(e),gmin(e));
                gmax(e) = max(g1(e),gmax(e));
                strainmin(e) = min(strain1(e),strainmin(e));
                strainmax(e) = max(strain1(e),strainmax(e));
              end
            end
          if(count==1)
             for i=1:ndof1
               ucmin(i,1) = ucvec(i,1);
               ucmax(i,1) = ucvec(i,1);
             end
           else %count>1
              for i=1:ndof1
                 if(ucvec(i,1)<ucmin(i,1))
                    ucmin(i,1) = ucvec(i,1);
                 end
                 if(ucvec(i,1)>ucmax(i,1))
                    ucmax(i,1) = ucvec(i,1);
                 end  
              end
          end
          if (mod(count,printfeq)==0)
        delt=toc/count;
        remaining=delt*(jmax-count);
        fprintf(1,"at %d estimated time remaining building data %s  delt %d\n", count,sec2hms(remaining),delt);
    end
end %end of force loop
end   %end of group loop
 
end   %end of element loop 
      uc = infsup(ucmin,ucmax);
%      for e=1:ne    % not sure why they are a matrix when only the only first index is 1
%       Cforce(:,e) = infsup(gmin(:,e),gmax(:,e));
%       Cstrain(:,e) = infsup(strainmin(:,e),strainmax(:,e)); %to obtain strains from comb.solution
%      end
     % replaced the above by the two lines below  RLM 9/9/22
      Cforce=infsup(gmin,gmax)  ;     
      Cstrain = infsup(strainmin,strainmax); 
   ComboTime=toc 
   fprintf(out,['vertex solution time  %s\n'],ComboTime);
 else
    fprintf(fid,'Combinatorial solution being skipped\n');
 end

fprintf(fid,'Nodal information - Assembled  model\n');
fprintf(fid,'Node   X      Y       Z       Restraints       Fx          Fy          Fz           beta-x    beta-y   beta-y\n');
for i=1:nnA
     fprintf(fid,'%2d   %4.1f   %4.1f   %4.1f    %d  %d  %d     %9.1f    %9.1f      %9.1f         %5.3f      %5.3f      %5.3f\n',i,xA(i),yA(i),zA(i),resxA(i),resyA(i),reszA(i),fxA(i),fyA(i),fzA(i),betaxA(i),betayA(i),betazA(i)); 
end
fprintf(fid,'Element information- Assembled model\n');
fprintf(fid,'Element    Nodes    Area            E        alfa      group no.    group gamma \n');
for i=1:ne
    fprintf(fid,'%2d        %2d  %2d    %f   %.3e   %5.3f         %2d           %5.3f\n',i,elementA(i,1),elementA(i,2),A(i),E(i),alfaA(i),gn(i), gammao(i));
end
fprintf(fid,'Groups Interval Solution\n');
fprintf(fid,'Number of Groups\n');
fprintf(fid,'%2d\n',ng);
fprintf(fid,'Group Number        Group Uncertainty   \n');
for i=1:ng
    fprintf(fid,'%2d                       %6.3f    \n',i,gamma1(i)); 
end
fprintf(fid,'Elemnt     Group \n');
mog = gn;
elng = find(mog);
mogg = nonzeros(mog);
lelg = length(elng);
for i=1:lelg
    fprintf(fid,'%2d          %2d          \n',elng(i),mogg(i)); 
end

fprintf(fid,'------------------------------------------------------------\n');
fprintf(fid,'\nNodal information - EBE model\n');
fprintf(fid,'Node   X      Y     Restraints       Fx          Fy           beta-x    beta-y\n');
for i=1:nn
     fprintf(fid,'%2d   %4.1f   %4.1f    %d    %d     %9.1f    %9.1f           %5.3f      %5.3f\n',i,xEBE(i),yEBE(i),resxEBE(i),resyEBE(i),fx(i),fy(i),betaxEBE(i),betayEBE(i)); 
end
fprintf(fid,'Element information- EBE model\n');
fprintf(fid,'Element    Nodes    Area            E        alfa      group no.    group gamma\n');
for i=1:ne
    fprintf(fid,'%2d        %2d  %2d    %f   %.3e   %5.3f         %2d           %5.3f\n',i,element(i,1),element(i,2),A(i),E(i),alfaA(i),gn(i), gammao(i));
end
fprintf(fid,'Groups Interval Solution\n');
fprintf(fid,'Number of Groups\n');
fprintf(fid,'%2d\n',ng);
fprintf(fid,'Group Number        Group Uncertainty   \n');
for i=1:ng
    fprintf(fid,'%2d                       %6.3f    \n',i,gamma1(i)); 
end
fprintf(fid,'Elemnt     Group \n');
mog = gn;
elng = find(mog);
mogg = nonzeros(mog);
lelg = length(elng);
for i=1:lelg
    fprintf(fid,'%2d          %2d          \n',elng(i),mogg(i)); 
end

fprintf(fid,'_______________________________________________________________________________________________________________________\n');
fprintf(fid,'Solution vector using EBE only\n');
fprintf(fid,'Element    node                     x- displacement                        y- displacement\n');
for e=1:ne
  n1 = 2*e-1; n2 = 2*e; %these are ebe nodes
  m1 = ebenode(n1); m2 = ebenode(n2); %these are the nodes in assembled model
  %changed output to file and not one  RLM 9/6/22
  fprintf(fid,'%2d          %d    [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,m1,inf(u(2*n1-1)),sup(u(2*n1-1)),inf(u(2*n1)),sup(u(2*n1)));
  fprintf(fid,'            %d    [%16.10e,%16.10e] [%16.10e,%16.10e]\n',m2,inf(u(2*n2-1)),sup(u(2*n2-1)),inf(u(2*n2)),sup(u(2*n2)));
end 
%----------------------------------------------------------------
fprintf(fid,'_______________________________________________________________________________________________________________________\n');
fprintf(fid,'Groups solution \n');
fprintf(fid,'Element    node                     x- displacement                        y- displacement\n');
for e=1:ne
  n1 = 2*e-1; n2 = 2*e; %these are ebe nodes
  m1 = ebenode(n1); m2 = ebenode(n2); %these are the nodes in assembled model
  fprintf(fid,'%2d          %d   [%16.10e,%16.10e] [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,m1,inf(us(3*n1-2)),sup(us(3*n1-2)),inf(us(3*n1-1)),sup(us(3*n1-1)),inf(us(2*n1)),sup(us(2*n1)));
  fprintf(fid,'            %d    [%16.10e,%16.10e] [%16.10e,%16.10e] [%16.10e,%16.10e]\n',m2,inf(us(3*n2-2)),sup(us(3*n2-2)),inf(us(3*n2-1)),sup(us(3*n2-1)),inf(us(3*n2)),sup(us(3*n2)));
end 
%----------------------------------------------------------------
if(option==1) %comb. solution is  chosen
    fprintf(fid,'_______________________________________________________________________________________________________________________\n');
    fprintf(fid,'\nSolution vector for all possible combinations\n');
    fprintf(fid,'Element    node                     x-displacement                            y-displacement\n');
      ucomb=infsup(zeros(1,2*nn),zeros(1,2*nn));  %rlm   set size of array for improving performance
      for e=1:ne
       n1 = 2*e-1;
       n2 = 2*e;  
       j1 = ebenode(n1);
       j2 = ebenode(n2);
        ucomb(3*n1-2) = uc(3*j1-2);ucomb(3*n1-1) = uc(3*j1-1); ucomb(3*n1)   = uc(3*j1); %at first node of element
       ucomb(3*n2-2) = uc(3*j2-2);ucomb(3*n2-1) = uc(3*j2-1); ucomb(3*n2)   = uc(3*j2); %at second node of element
       %changed output to file and not one  RLM 9/6/22
       fprintf(fid,'%2d          %d    [%16.10e,%16.10e] [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,j1,inf(ucomb(3*n1-2)),sup(ucomb(3*n1-2)),inf(ucomb(3*n1-1)),sup(ucomb(3*n1-1)),inf(ucomb(3*n1)),sup(ucomb(3*n1)));
       fprintf(fid,'            %d    [%16.10e,%16.10e] [%16.10e,%16.10e] [%16.10e,%16.10e]\n',j2,inf(ucomb(3*n2-2)),sup(ucomb(3*n2-2)),inf(ucomb(3*n2-1)),sup(ucomb(3*n2-1)),inf(ucomb(3*n2)),sup(ucomb(2*n2)));
     end
end
if(option==2)
  fprintf(fid,'Combinatorial solution is skipped\n');
end
if(option==1)
  fprintf(fid,'_______________________________________________________________________________________________________________________\n');
  fprintf(fid,'WIDTH AND ERROR OF SOLUTION VECTOR\n');

  fprintf(fid,'-----------------------------------------------------------------------------------------------------------------------\n');
  fprintf(fid,'             Combin.solution      New approach       Error%%       LB-Error%%    UB-Error%%\n');
  fprintf(fid,'________________________________________________________________________________________________________________________\n');
  fprintf(fid,'Elem node                    \n');
for e=1:ne
    n1  = 2*e-1; n2  = 2*e ;
    m1 = ebenode(n1); m2 = ebenode(n2); 
    d1x = 2*n1-1; %dof1
    d1y = 2*n1;   %dof2
    %d1z = 3*n1
    d2x = 2*n2-1; %dof3
    d2y = 2*n2;   %dof4
    %d2z = 3*n2;
    %widths for present approach  (Change to group results RLM 9/8/22)
    numerb1 = sup(us(d1x))-inf(us(d1x));%dof1
    numerb2 = sup(us(d1y))-inf(us(d1y));%dof2
    numerb3 = sup(us(d2x))-inf(us(d2x));%dof3
    numerb4 = sup(us(d2y))-inf(us(d2y));%dof4
    if(option==1) %comb. solution is chosen
      %widths for combinatorial solution
      denom1 = sup(ucomb(d1x))-inf(ucomb(d1x)); %dof1
      denom2 = sup(ucomb(d1y))-inf(ucomb(d1y)); %dof2
      denom3 = sup(ucomb(d2x))-inf(ucomb(d2x)); %dof3
      denom4 = sup(ucomb(d2y))-inf(ucomb(d2y)); %dof4
    else 
      denom1 = 0;
      denom2 = 0;
      denom3 = 0;
      denom4 = 0;
    end
    if (denom1>1.0e-15&&numerb1>1.0e-15)
      error2=100*numerb1/denom1-100; %for present approach
      diff1 = 100*abs((inf(ucomb(d1x))-inf(u(d1x)))/inf(ucomb(d1x))); 
      diff2 = 100*abs((sup(ucomb(d1x))-sup(u(d1x)))/sup(ucomb(d1x))); 
      fprintf(fid,'%2d    %d  u  %17.10e   %17.10e  %9.6f     %7.4f      %7.4f\n',e,m1,denom1,numerb1,error2,diff1,diff2);
    else
      fprintf(fid,'%2d    %d  u  %17.10e   %17.10e      ----        ----         ----\n',e,m1,denom1,numerb1);
    end
    if (denom2>1.0e-15&&numerb2>1.0e-15)
      error2=100*numerb2/denom2-100; %for present approach
      diff1 = 100*abs((inf(ucomb(d1y))-inf(us(d1y)))/inf(ucomb(d1y))); 
      diff2 = 100*abs((sup(ucomb(d1y))-sup(us(d1y)))/sup(ucomb(d1y))); 
      fprintf(fid,'         v  %17.10e   %17.10e  %9.6f     %7.4f      %7.4f\n',denom2,numerb2,error2,diff1,diff2);
    else
      fprintf(fid,'         v  %17.10e   %17.10e      ----        ----         ----\n',denom2,numerb2);
    end
    if (denom3>1.0e-15&&numerb3>1.0e-15)
      error2=100*numerb3/denom3-100; %for present approach
      diff1 = 100*abs((inf(ucomb(d2x))-inf(us(d2x)))/inf(ucomb(d2x))); 
      diff2 = 100*abs((sup(ucomb(d2x))-sup(us(d2x)))/sup(ucomb(d2x))); 
      fprintf(fid,'      %d  u  %17.10e   %17.10e  %9.6f     %7.4f      %7.4f\n',m2,denom3,numerb3,error2,diff1,diff2);
    else
      fprintf(fid,'      %d  u  %17.10e   %17.10e     ----         ----         ----\n',m2,denom3,numerb3);
    end
    if (denom4>1.0e-15&&numerb4>1.0e-15)
      error2=100*numerb4/denom4-100; %for present approach
      diff1 = 100*abs((inf(ucomb(d2y))-inf(us(d2y)))/inf(ucomb(d2y))); 
      diff2 = 100*abs((sup(ucomb(d2y))-sup(us(d2y)))/sup(ucomb(d2y))); 
      fprintf(fid,'         v  %17.10e   %17.10e  %9.6f     %7.4f      %7.4f\n',denom4,numerb4,error2,diff1,diff2);
    else
      fprintf(fid,'         v  %17.10e   %17.10e     ----         ----         ----\n',denom4,numerb4);
    end
    fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
end 
   
end
fprintf(fid,'_________________________________________________________________________________________________\n');
fprintf(fid,'A X I A L  F O R C E   I N   M E M B E R S \n');
fprintf(fid,'---------------------------------------------------------------------------------------------------\n');
fprintf(fid,'MEMBER  Combinatorial Solution      Present solution       Error%%       LB-Error%%     UB-Error%%\n');
for i=1:ne
  w1 = sup(eforces(i,1))-inf(eforces(i,1));
  if(option==1)
     w2 = sup(Cforce(i))-inf(Cforce(i)); %comb solution is chosen
     if(inf(Cforce(i))>0)
       diff1 = 100*(inf(Cforce(i))-inf(eforces(i,1)))/inf(Cforce(i));
     end
     if(inf(Cforce(i))<0)
       diff1 = -100*(inf(Cforce(i))-inf(eforces(i,1)))/inf(Cforce(i));
     end
     if(sup(Cforce(i))>0)
       diff2 = 100*(sup(eforces(i,1))-sup(Cforce(i)))/sup(Cforce(i));
     end
     if(sup(Cforce(i))<0)
       diff2 = -100*(sup(eforces(i,1))-sup(Cforce(i)))/sup(Cforce(i));
     end
     if(w1>0&&w2>0)
         error = 100*(w1/w2)-100.0;
          fprintf(fid,'%3d  [%10.6f,%10.6f]  [%10.6f,%10.6f]   %10.6f    %10.6f    %10.6f\n',i,inf(Cforce(i)),sup(Cforce(i)),inf(eforce(i,1)),sup(eforce(i,1)),error,diff1,diff2);
     else
         fprintf(fid,'%3d  [%10.6f,%10.6f]   [%10.6f,%10.6f]      ----         ----        -----\n',i,inf(Cforce(i)),sup(Cforce(i)),inf(eforce(i,1)),sup(eforce(i,1)));
     end
  else
      w2 = 0.0;  %comb solution is skipped
      error = 0.0;
      fprintf(fid,'%3d       -----              [%15.6f,%15.6f]         ----         ----        -----\n',i,inf(eforce(i,1)),sup(eforce(i,1)));
  end
end %of for loop on elements
%------------------------------------------------------------
%Display group forces
%------------------------------------------------------------
fprintf(fid,'_________________________________________________________________________________________________\n');
fprintf(fid,'G R O U P S  A X I A L  F O R C E   I N   M E M B E R S \n');
fprintf(fid,'---------------------------------------------------------------------------------------------------\n');
fprintf(fid,'MEMBER  Combinatorial Solution      Present solution       Error%%       LB-Error%%     UB-Error%%\n');
for i=1:ne
  w1 = sup(eforces(i,1))-inf(eforces(i,1));
  if(option==1)
     w2 = sup(Cforce(i))-inf(Cforce(i)); %comb solution is chosen
     if(inf(Cforce(i))>0)
       diff1 = 100*(inf(Cforce(i))-inf(eforces(i,1)))/inf(Cforce(i));
     end
     if(inf(Cforce(i))<0)
       diff1 = -100*(inf(Cforce(i))-inf(eforces(i,1)))/inf(Cforce(i));
     end
     if(sup(Cforce(i))>0)
       diff2 = 100*(sup(eforces(i,1))-sup(Cforce(i)))/sup(Cforce(i));
     end
     if(sup(Cforce(i))<0)
       diff2 = -100*(sup(eforces(i,1))-sup(Cforce(i)))/sup(Cforce(i));
     end
     if(w1>0&&w2>0)
         error = 100*(w1/w2)-100.0;
          fprintf(fid,'%3d  [%10.6f,%10.6f]  [%10.6f,%10.6f]   %10.6f    %10.6f    %10.6f\n',i,inf(Cforce(i)),sup(Cforce(i)),inf(eforces(i,1)),sup(eforces(i,1)),error,diff1,diff2);
     else
         fprintf(fid,'%3d  [%10.6f,%10.6f]   [%10.6f,%10.6f]      ----         ----        -----\n',i,inf(Cforce(i)),sup(Cforce(i)),inf(eforces(i,1)),sup(eforces(i,1)));
     end
  else
      w2 = 0.0;  %comb solution is skipped
      error = 0.0;
      fprintf(fid,'%3d       -----              [%15.6f,%15.6f]         ----         ----        -----\n',i,inf(eforces(i,1)),sup(eforces(i,1)));
  end
end %of for loop on elements


%------------------------------------------------------------
%changed output to file and not one  RLM 9/6/22
fprintf(fid,'_____________________________________________________________________________\n');
fprintf(fid,'S T R E S S E S   A N D   S T R A I N S    I N   M E M B E R S \n');
fprintf(fid,'-----------------------------------------------------------------------------\n');
fprintf(fid,'MEMBER           stress (N/m^2)                 Strain (computed from stress) \n');
for e=1:ne    
    estress(e) = eforces(e,1)/(A(e)*Gamma(e));  % RLM 9/9/22 changed eforce to eforces to include groups
    estrain(e) = estress(e)/(Alfa(e)*E(e)) ;  % RLM 9/6/22  added Group term for strain
    %changed output to file and not one  RLM 9/6/22
     fprintf(fid,'%2d      [%17.8e,%17.8e]   [%17.8e,%17.8e]\n',e,inf(estress(e)),sup(estress(e)),inf(estrain(e)),sup(estrain(e)));
end
%extracting stresses in elements
for i=ndofnew-ne+1:ndofnew
    j=i-ndof-ndof1-nlamda-ne;
    strain(j) = u1s(i); % changed u1 to u1s  RLM 9/8/22
end
 %changed output to file and not one  RLM 9/6/22

fprintf(fid,'___________________________________________________________________________________________________\n');
fprintf(fid,'S T R A I N S   I N   M E M B E R S \n');
fprintf(fid,'---------------------------------------------------------------------------------------------------\n');
fprintf(fid,'MEMBER  Combinatorial Solution              New Formulation                   Error%%      LB-Error%%     UB-Error%%\n');
for e=1:ne
  
  
  if(option==1)
      %moved to below if statement so not report if no Combinaorial solution RLM 9/7/22  
      w1 = sup(Cstrain(e))-inf(Cstrain(e));
     w2 = sup(strain(e))-inf(strain(e)); %comb solution is chosen
     if(inf(Cstrain(e))>0)
       diff1 = 100*(inf(Cstrain(e))-inf(strain(e)))/inf(Cstrain(e));
     end
     if(inf(Cstrain(e))<0)
       diff1 = -100*(inf(Cstrain(e))-inf(strain(e)))/inf(Cstrain(e));
     end
     if(sup(Cstrain(e))>0)
       diff2 = 100*(sup(strain(e))-sup(Cstrain(e)))/sup(Cstrain(e));
     end
     if(sup(Cstrain(e))<0)
       diff2 = -100*(sup(strain(e))-sup(Cstrain(e)))/sup(Cstrain(e));
     end
     
     
     if(w2>0)  % RLM maybe check should be w1 since combinations
          error = 100*(w2/w1)-100.0;
          fprintf(fid,'%3d  [%15.6e,%15.6e]  [%15.6e,%15.6e]   %10.6f    %10.6f    %10.6f\n',e,inf(Cstrain(e)),sup(Cstrain(e)),inf(strain(e)),sup(strain(e)),error,diff1,diff2);
     else
         fprintf(fid,'%3d  [%15.6e,%15.6e]   [%15.6e,%15.6e]      ----         ----        -----\n',e,inf(Cstrain(e)),sup(Cstrain(e)),inf(strain(e)),sup(strain(e)));
     end
  else
      w2 = 0.0;  %comb solution is skipped
      error = 0.0;
      fprintf(fid,'%3d       -----              [%15.6e,%15.6e]         ----         ----        -----\n',e,inf(strain(e)),sup(strain(e)));
  end
end %of for loop on elements
toc

