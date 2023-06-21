function [option,uc,Cforce,Cstrain ]=groupcomb3 (fid,ne,ng,nnA,nforce,E,A,Alfa,Gamma,ForceA,xA,yA,zA,elementA,gn,resxA,resyA,reszA);
% program to solve the combination solution for a 3 d linear truss
ndof1=3*nnA;
for i=1:nnA
      if(resxA(i)==0) %x-dof is restrained
        bnA = bnA+1;
        ifixA(bnA) = 3*i-2;
    end
   if(resyA(i)==0) %y-dof is restrained
        bnA = bnA+1;
        ifixA(bnA) = 3*i-1;
   end
   if(reszA(i)==0) %y-dof is restrained
        bnA = bnA+1;
        ifixA(bnA) = 3*i;
   end
end
j=0;
 intforce = 0;
 for i=1:ndof1
    if(mid(ForceA(i)~=0))
        j=j+1;
        loc(j)=i;
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
     % bin1=ones(intforce,1);
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
             % buf = int2str(bin1)';
     % ncombF= bintodec(buf);
     % bin2 = ones(nE,1);
     % buf = int2str(bin2)';
     % ncombE= bintodec(buf);
     % %added bin3 for group data
     % bin3 = ones(ng1,1);
     % buf = int2str(bin3)';
     % ncombG= bintodec(buf);
     fprintf(fid,'Number of combinations force %d, combinations E %d, combinations A (groups) %d\n',ncombF+1,ncombE+1,ncombG+1);
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
                j=loc(i); %location of dof carrying non-zero load
               if(combF(i)==0)  %set load to lower bound of ForceA
                    P1(j)= inf(ForceA(j));
               else
                    P1(j)=sup(ForceA(j));%set load to upper bound of ForceA
               end
            end
            Kic=zeros(ndof1,ndof1);
            for e=1:ne
             conn1=elementA(e,:);
             dof2=[3*conn1(1)-2 3*conn1(1)-1 3*conn1(1) 3*conn1(2)-2 3*conn1(2)-1 3*conn1(2)];
             x1=xA(conn1(1));
             y1=yA(conn1(1));
             z1=zA(conn1(1));
             x2=xA(conn1(2));
             y2=yA(conn1(2));
             z2=zA(conn1(2));
             le=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
             cx=(x2-x1)/le;
             cy=(y2-y1)/le; 
             cz=(z2-z1)/le;
             % change to A1 for groups
             kc=A1(e)*(E1(e))/le;
             ke=kc*[cx*cx  cx*cy  cx*cz -cx*cx -cx*cy  -cx*cz;
                    cy*cx  cy*cy  cy*cz -cy*cx -cy*cy  -cy*cz;      
                    cz*cx  cz*cy  cz*cz -cz*cx -cz*cy  -cz*cz;
                    -cx*cx  -cx*cy  -cx*cz cx*cx cx*cy  cx*cz;
                    -cy*cx  -cy*cy  -cy*cz cy*cx cy*cy  cy*cz;      
                    -cz*cx  -cz*cy  -cz*cz cz*cx cz*cy  cz*cz];
                for i1=1:6
                j1=dof2(i1);
                for i2=1:66
                    j2=dof2(i2);
            Kic(j1,j2)=Kic(j1,j2)+ke(i1,i2);        
                end
            end
            end
            Kic(ifixA,:)=0;                   
            Kic(:,ifixA)=0;                    
            Kic(ifixA,ifixA)=eye(length(ifixA));
            
            ucvec =Kic\P1';
 
            for e=1:ne
              
             conn1=elementA(e,:);
             dof2=[3*conn1(1)-2 3*conn1(1)-1 3*conn1(1) 3*conn1(2)-2  3*conn1(2)-1 3*conn1(2)];
             x1=xA(conn1(1));
             y1=yA(conn1(1));
             z1=zA(conn1(1));
             x2=xA(conn1(2));
             y2=yA(conn1(2));
             z2=zA(conn1(2));
             le=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
             cx=(x2-x1)/le;
             cy=(y2-y1)/le; 
             cz=(z2-z1)/le;
                        
             eps=[-cx -cy -cz cx cy cz]/le*ucvec(dof2);
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
   fprintf(fid,['vertex solution time  %s\n'],ComboTime);
 else
    fprintf(fid,'Combinatorial solution being skipped\n');
 end


