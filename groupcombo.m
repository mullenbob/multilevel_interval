function [option,uc,Cforce,Cstrain ]=groupcombo (fid,ne,ng,nnA,nforce,E,A,Alfa,Gamma,gamma,ForceA,xA,yA,elementA,gn,resxA,resyA)
% program to solve the combination solution for a 2 d linear truss
% return values are option =1 solution option = 2 no solution
% uc is interval node dis[lacemend
% Cforce is interval force
% Cstrain is interval strain
% input values
% fid           file for display
% ne            number of elements
% ng            number of groups
% nnA           number of nodes in assembled model
% nforce
% E             element midpoint modulus
% A             element midpoint area
% Alfa          element interval multipler
% Gamma         interval multiplier for each element
% gamma         interval radius for each group
% ForceA        interval Applied force to assembled model
% xA            x coordinate for each node
% yA            y coordinate  for each node
% elementA      global node numbers for each element    
%  gn           group number for each element
% resxA         x direction bc code  = 0 fixed =1 free
% resyA         y direction bc code 
% calculated number of assembled dof
ndof1=2*nnA;
% construct an index array for fixed nodes   local variable ifixA
bnA=0;
for i=1:nnA
    if(resxA(i)==0) %x-dof is restrained
        bnA = bnA+1;
        ifixA(bnA) = 2*i-1;
    end
   if(resyA(i)==0) %y-dof is restrained
        bnA = bnA+1;
        ifixA(bnA) = 2*i;
   end
end
% count the number of independent interval  forces

 intforce = 0;
 for i=1:ndof1
    if((sup(ForceA(i))-inf(ForceA(i)))>0)
        intforce = intforce+1;
        loc(intforce)=i;    
    end
 end
 ng1 = 0;
 for i = 1:ng
     if (gamma(i))>0
         ng1=ng1+1;
     end
 end
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
     
     ncombF=2^intforce;
     ncombE=2^ne;
     ncombG=2^ng1;
     % allocate space for results
             gmin=zeros(ne,1);
             gmax=gmin;
             strainmin=gmin;
             strainmax=gmin;
             ucmin=zeros(ndof1,1);
             ucmax=ucmin;
     
     fprintf(fid,'Number of combinations force %d, combinations E %d, combinations G (groups) %d\n',ncombF,ncombE,ncombG);
     if (ng1==0)
     ncombG=0;
     end

     if(ne==0)
         ncombE=0;
     end
     if(intforce==0)
         ncombF=0;
     end
     %generating combinations of upper and lower bounds of loads and E and
     %and groups 
     count=0;
E1=E;
for ncE=0:ncombE-1
         combE =dectobin(ncE,ne);
         for i=1:ne
            if(combE(i)==0)
                E1(i)=Q1(i,1);
            else
                E1(i)=Q1(i,2);
            end
         end

for ncG=0:ncombG-1
             combG =dectobin(ncG,ng)';
             %modify area for group uncertainty RLM 9/8/22
             A1=A;
         for i=1:ne
             ii=gn(i);
             
             if (ii > 0)
            if(combG(ii)==0)
                A1(i)=A(i)*Q2(i,1);
            else
                A1(i)=A(i)*Q2(i,2);
            end
             end
         end


for ncF=0:ncombF-1
            count = count+1;
            combF =dectobin(ncF,intforce)';
            % set forces to mid point and correct for interval values latter
            % for i=1:ndof1
                P1=mid(ForceA);
            
            for i=1:intforce
                j=loc(i); %location of dof carrying non-interval non-zero  load
               if(combF(i)==0)  %set load to lower bound of ForceA
                    P1(j)= inf(ForceA(j));
               else
                    P1(j)=sup(ForceA(j));%set load to upper bound of ForceA
               end
            end
            Kic=zeros(ndof1,ndof1);
            for e=1:ne
             conn1=elementA(e,:);
             dof2=[2*conn1(1)-1 2*conn1(1) 2*conn1(2)-1 2*conn1(2)];
             x1=xA(conn1(1));
             y1=yA(conn1(1));
             x2=xA(conn1(2));
             y2=yA(conn1(2));
             le=sqrt((x2-x1)^2+(y2-y1)^2);
             c=(x2-x1)/le;
             s=(y2-y1)/le; 
            
             % change to A1 for groups
             kc=A1(e)*(E1(e))/le;
             ke=kc*[c*c  c*s  -c*c -c*s;
                    c*s  s*s  -c*s -s*s;
                    -c*c -c*s  c*c  c*s;
                    -c*s -s*s  c*s  s*s];
            for i1=1:4
                j1=dof2(i1);
                for i2=1:4
                    j2=dof2(i2);
            Kic(j1,j2)=Kic(j1,j2)+ke(i1,i2);        
                end
            end
            end
            Kic(ifixA,:)=0;                   
            Kic(:,ifixA)=0;                    
            Kic(ifixA,ifixA)=eye(length(ifixA));
            
            ucvec =Kic\P1;
 
            for e=1:ne
              
             conn1=elementA(e,:);
             dof2=[2*conn1(1)-1 2*conn1(1) 2*conn1(2)-1 2*conn1(2)];
             x1=xA(conn1(1));
             y1=yA(conn1(1));
             x2=xA(conn1(2));
             y2=yA(conn1(2));
             le=sqrt((x2-x1)^2+(y2-y1)^2);
             c=(x2-x1)/le;
             s=(y2-y1)/le; 
             
             eps=[-c -s c s]/le*ucvec(dof2);
             %changed to use A1 to account for groups
              
              g1=eps*A1(e)*E1(e);
              if(count==1)
                  gmin(e)= g1;
                  gmax(e)= g1;
                  strainmin(e) = eps;
                  strainmax(e) = eps;
              else
                gmin(e) = min(g1,gmin(e));
                gmax(e) = max(g1,gmax(e));
                strainmin(e) = min(eps,strainmin(e));
                strainmax(e) = max(eps,strainmax(e));
              end
            end
          if(count==1)
             ucmin=ucvec;
             ucmax=ucvec;
           else %count>1
              for i=1:ndof1
                 if(ucvec(i)<ucmin(i))
                    ucmin(i) = ucvec(i);
                 end
                 if(ucvec(i)>ucmax(i))
                    ucmax(i) = ucvec(i);
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
   ComboTime=toc ;
   fprintf(fid,'vertex solution time  %s\n',ComboTime);
 else
    fprintf(fid,'Combinatorial solution being skipped\n');
    %set dummy return values
    uc=0;
    Cforce=0;
    Cstrain=0;
 end


