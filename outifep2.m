  
function outifep2(fid,strainflag,option,nnA,ndof,ndof1,ndofnew,ng,gamma1,gn,ne,elementA,Alfa,E,A,Gamma,uu,u1,us,ugs,ucombo,Cforce,Cstrain)
% fprintf(fid,'------------------------------------------------------------\n');
% fprintf(fid,'\nNodal information - EBE model\n');
% fprintf(fid,'Node   X      Y     Restraints       Fx          Fy           beta-x    beta-y\n');
% for i=1:nn
%      fprintf(fid,'%2d   %4.1f   %4.1f    %d    %d     %9.1f    %9.1f           %5.3f      %5.3f\n',i,xEBE(i),yEBE(i),resxEBE(i),resyEBE(i),fx(i),fy(i),betaxEBE(i),betayEBE(i)); 
% end
% fprintf(fid,'Element information- EBE model\n');
% fprintf(fid,'Element    Nodes    Area            E        alfa      group no.    group gamma\n');
% for i=1:ne
%     fprintf(fid,'%2d        %2d  %2d    %f   %.3e   %5.3f         %2d           %5.3f\n',i,element(i,1),element(i,2),A(i),E(i),alfaA(i),gn(i), gammao(i));
% end

% print out the centered solution
fprintf(fid,'_______________________________________________________________________________________________________________________\n');
fprintf(fid,'Centered Solution \n');
fprintf(fid,'Element    node x- displacement                 y- displacement\n');
for e=1:ne
  
  m1 = elementA(e,1); m2=elementA(e,2);
  
  fprintf(fid,'%2d          %d    %16.10e             %16.10e\n',e,m1,uu(ndof+2*m1-1),uu(ndof+2*m1));
  fprintf(fid,'%2d          %d    %16.10e             %16.10e \n',e,m2,uu(ndof+2*m2-1),uu(ndof+2*m2));
  
end 
fprintf(fid,'_______________________________________________________________________________________________________________________\n');
fprintf(fid,'Solution vector using EBE old method\n');
fprintf(fid,'Element    node                     x- displacement                        y- displacement\n');
for e=1:ne
  m1 = elementA(e,1); m2=elementA(e,2);
  fprintf(fid,'%2d          %d    [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,m1,inf(u1(ndof+2*m1-1)),sup(u1(ndof+2*m1-1)),inf(u1(ndof+2*m1)),sup(u1(ndof+2*m1)));
  fprintf(fid,'%2d          %d    [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,m2,inf(u1(ndof+2*m2-1)),sup(u1(ndof+2*m2-1)),inf(u1(ndof+2*m2)),sup(u1(ndof+2*m2)));
end 
fprintf(fid,'_______________________________________________________________________________________________________________________\n');
fprintf(fid,'Solution vector using EBE only union method\n');
fprintf(fid,'Element    node                     x- displacement                        y- displacement\n');
for e=1:ne
  m1 = elementA(e,1); m2=elementA(e,2);
  fprintf(fid,'%2d          %d    [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,m1,inf(ugs(ndof+2*m1-1)),sup(ugs(ndof+2*m1-1)),inf(ugs(ndof+2*m1)),sup(ugs(ndof+2*m1)));
  fprintf(fid,'%2d          %d    [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,m2,inf(ugs(ndof+2*m2-1)),sup(ugs(ndof+2*m2-1)),inf(ugs(ndof+2*m2)),sup(ugs(ndof+2*m2)));
end 
%----------------------------------------------------------------
fprintf(fid,'_______________________________________________________________________________________________________________________\n');
fprintf(fid,'Groups solution \n');
fprintf(fid,'Element    node                     x- displacement                        y- displacement\n');
for e=1:ne
  m1 = elementA(e,1); m2=elementA(e,2);
  fprintf(fid,'%2d          %d    [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,m1,inf(us(ndof+2*m1-1)),sup(us(ndof+2*m1-1)),inf(us(ndof+2*m1)),sup(us(ndof+2*m1)));
  fprintf(fid,'%2d          %d    [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e,m2,inf(us(ndof+2*m2-1)),sup(us(ndof+2*m2-1)),inf(us(ndof+2*m2)),sup(us(ndof+2*m2)));
end
%----------------------------------------------------------------
if (option == 1)
fprintf(fid,'_______________________________________________________________________________________________________________________\n');
fprintf(fid,'end-point solution  \n');
fprintf(fid,'    node                     x- displacement                        y- displacement\n');
for e=1:nnA 
  fprintf(fid,'%2d     [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e, inf(ucombo(2*e-1)),sup(ucombo(2*e-1)),inf(ucombo(2*e)),sup(ucombo(2*e)));
end 
fprintf(fid,'\n end-point solution  \n');
fprintf(fid,'   element                   force                     strain \n');
for e=1:ne
  fprintf(fid,'%2d      [%16.10e,%16.10e] [%16.10e,%16.10e]\n',e, inf(Cforce(e)),sup(Cforce(e)),inf(Cstrain(e)),sup(Cstrain(e)));
end 
end
  fprintf(fid,'_______________________________________________________________________________________________________________________\n');
  fprintf(fid,'WIDTH AND ERROR OF SOLUTION VECTOR\n');

  fprintf(fid,'-----------------------------------------------------------------------------------------------------------------------\n');
  fprintf(fid,'             Combin.solution      New approach       Error%%       LB-Error%%    UB-Error%%\n');
  fprintf(fid,'________________________________________________________________________________________________________________________\n');
  fprintf(fid,'Elem node                    \n');
for e=1:ne
   
    m1 = elementA(e,1); m2=elementA(e,2);
    
%     d1x = 2*n1-1; %dof1
%     d1y = 2*n1;   %dof2
%     d2x = 2*n2-1; %dof3
%     d2y = 2*n2;   %dof4
    d1x = 2*ne+2*m1-1; %dof1
    d1y = 2*ne+2*m1;   %dof2
    d2x = 2*ne+2*m2-1; %dof3
    d2y = 2*ne+2*m2;   %dof4
    %widths for present approach  (Change to group results RLM 9/8/22)
    numerb1 = sup(us(d1x))-inf(us(d1x));%dof1
    numerb2 = sup(us(d1y))-inf(us(d1y));%dof2
    numerb3 = sup(us(d2x))-inf(us(d2x));%dof3
    numerb4 = sup(us(d2y))-inf(us(d2y));%dof4
    if(option==1) %comb. solution is chosen
      %widths for combinatorial solution
      denom1 = sup(ucombo(d1x-ndof))-inf(ucombo(d1x-ndof)); %dof1
      denom2 = sup(ucombo(d1y-ndof))-inf(ucombo(d1y-ndof)); %dof2
      denom3 = sup(ucombo(d2x-ndof))-inf(ucombo(d2x-ndof)); %dof3
      denom4 = sup(ucombo(d2y-ndof))-inf(ucombo(d2y-ndof)); %dof4
    else 
      denom1 = 0;
      denom2 = 0;
      denom3 = 0;
      denom4 = 0;
    end
      
    if (denom1>1.0e-15&&numerb1>1.0e-15)
      error2=100*numerb1/denom1-100; %for present approach
      diff1 = 100*abs((inf(ucombo(d1x-ne*2))-inf(us(d1x)))/inf(ucombo(d1x-ne*2))); 
      diff2 = 100*abs((sup(ucombo(d1x-ne*2))-sup(us(d1x)))/sup(ucombo(d1x-ne*2))); 
      fprintf(fid,'%2d    %d  u  %17.10e   %17.10e  %9.6f     %7.4f      %7.4f\n',e,m1,denom1,numerb1,error2,diff1,diff2);
    else
      fprintf(fid,'%2d    %d  u  %17.10e   %17.10e      ----        ----         ----\n',e,m1,denom1,numerb1);
    end
    if (denom2>1.0e-15&&numerb2>1.0e-15)
      error2=100*numerb2/denom2-100; %for present approach
      diff1 = 100*abs((inf(ucombo(d1y-ne*2))-inf(us(d1y)))/inf(ucombo(d1y-ne*2))); 
      diff2 = 100*abs((sup(ucombo(d1y-ne*2))-sup(us(d1y)))/sup(ucombo(d1y-ne*2))); 
      fprintf(fid,'         v  %17.10e   %17.10e  %9.6f     %7.4f      %7.4f\n',denom2,numerb2,error2,diff1,diff2);
    else
      fprintf(fid,'         v  %17.10e   %17.10e      ----        ----         ----\n',denom2,numerb2);
    end
    if (denom3>1.0e-15&&numerb3>1.0e-15)
      error2=100*numerb3/denom3-100; %for present approach
      diff1 = 100*abs((inf(ucombo(d2x-ne*2))-inf(us(d2x)))/inf(ucombo(d2x-ne*2))); 
      diff2 = 100*abs((sup(ucombo(d2x-ne*2))-sup(us(d2x)))/sup(ucombo(d2x-ne*2))); 
      fprintf(fid,'      %d  u  %17.10e   %17.10e  %9.6f     %7.4f      %7.4f\n',m2,denom3,numerb3,error2,diff1,diff2);
    else
      fprintf(fid,'      %d  u  %17.10e   %17.10e     ----         ----         ----\n',m2,denom3,numerb3);
    end
    if (denom4>1.0e-15&&numerb4>1.0e-15)
      error2=100*numerb4/denom4-100; %for present approach
      diff1 = 100*abs((inf(ucombo(d2y-ne*2))-inf(us(d2y)))/inf(ucombo(d2y-ne*2))); 
      diff2 = 100*abs((sup(ucombo(d2y-ne*2))-sup(us(d2y)))/sup(ucombo(d2y-ne*2))); 
      fprintf(fid,'         v  %17.10e   %17.10e  %9.6f     %7.4f      %7.4f\n',denom4,numerb4,error2,diff1,diff2);
    else
      fprintf(fid,'         v  %17.10e   %17.10e     ----         ----         ----\n',denom4,numerb4);
    end
    fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
end 
   

fprintf(fid,'_________________________________________________________________________________________________\n');
fprintf(fid,'A X I A L  F O R C E   I N   M E M B E R S  EBE only \n');
fprintf(fid,'---------------------------------------------------------------------------------------------------\n');
fprintf(fid,'MEMBER  Combinatorial Solution      Present solution       Error%%       LB-Error%%     UB-Error%%\n');
for i=1:ne
    i2=ndof+ndof1+2*i-1;
  w1 = sup(u1(i2))-inf(u1(i2));
  if(option==1)
     w2 = sup(Cforce(i))-inf(Cforce(i)); %comb solution is chosen
     if(inf(Cforce(i))>0)
       diff1 = 100*(inf(Cforce(i))-inf(u1(i2)))/inf(Cforce(i));
     end
     if(inf(Cforce(i))<0)
       diff1 = -100*(inf(Cforce(i))-inf(u1(i2)))/inf(Cforce(i));
     end
     if(sup(Cforce(i))>0)
       diff2 = 100*(sup(u1(i2))-sup(Cforce(i)))/sup(Cforce(i));
     end
     if(sup(Cforce(i))<0)
       diff2 = -100*(sup(u1(i2))-sup(Cforce(i)))/sup(Cforce(i));
     end
     
     if(w1>0&&w2>0)
         error = 100*(w1/w2)-100.0;
          fprintf(fid,'%3d  [%10.6f,%10.6f]  [%10.6f,%10.6f]   %10.6f    %10.6f    %10.6f\n',i,inf(Cforce(i)),sup(Cforce(i)),inf(u1(i2)),sup(u1(i2)),error,diff1,diff2);
     else
         fprintf(fid,'%3d  [%10.6f,%10.6f]   [%10.6f,%10.6f]      ----         ----        -----\n',i,inf(Cforce(i)),sup(Cforce(i)),inf(u1(i2)),sup(u1(i2)));
     end
  else
      
      fprintf(fid,'%3d       -----              [%15.6f,%15.6f]         ----         ----        -----\n',i,inf(us(ndof+ndof1+i)),sup(us(ndof+ndof1+i)));
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
    i2=ndof+ndof1+2*i-1;
  w1 = sup(us(i2))-inf(us(i2));
  if(option==1)
     w2 = sup(Cforce(i))-inf(Cforce(i)); %comb solution is chosen
     if(inf(Cforce(i))>0)
       diff1 = 100*(inf(Cforce(i))-inf(us(i2)))/inf(Cforce(i));
     end
     if(inf(Cforce(i))<0)
       diff1 = -100*(inf(Cforce(i))-inf(us(i2)))/inf(Cforce(i));
     end
     if(sup(Cforce(i))>0)
       diff2 = 100*(sup(us(i2))-sup(Cforce(i)))/sup(Cforce(i));
     end
     if(sup(Cforce(i))<0)
       diff2 = -100*(sup(us(i2))-sup(Cforce(i)))/sup(Cforce(i));
     end
     if(w1>0&&w2>0)
         error = 100*(w1/w2)-100.0;
          fprintf(fid,'%3d  [%10.6f,%10.6f]  [%10.6f,%10.6f]   %10.6f    %10.6f    %10.6f\n',i,inf(Cforce(i)),sup(Cforce(i)),inf(us(i2)),sup(us(i2)),error,diff1,diff2);
     else
         fprintf(fid,'%3d  [%10.6f,%10.6f]   [%10.6f,%10.6f]      ----         ----        -----\n',i,inf(Cforce(i)),sup(Cforce(i)),inf(us(i2)),sup(us(i2)));
     end
  else
      
      fprintf(fid,'%3d       -----              [%15.6f,%15.6f]         ----         ----        -----\n',i,inf(us(ndof+ndof1+i)),sup(us(ndof+ndof1+i)));
  end
end %of for loop on elements


%------------------------------------------------------------
%changed output to file and not one  RLM 9/6/22
fprintf(fid,'_____________________________________________________________________________\n');
fprintf(fid,'S T R E S S E S   A N D   S T R A I N S    I N   M E M B E R S \n');
fprintf(fid,'-----------------------------------------------------------------------------\n');
fprintf(fid,'MEMBER           stress (N/m^2)                 Strain (computed from stress) \n');
estress=zeros(ne,1)*infsup(1.,1.);
estrain=estress;
for e=1:ne    
    estress(e) = us(ndof+ndof1+2*e-1)/(A(e)*Gamma(e));  % RLM 9/9/22 changed eforce to eforces to include groups
    estrain(e) = estress(e)/(Alfa(e)*E(e)) ;  % RLM 9/6/22  added Group term for strain
    %changed output to file and not one  RLM 9/6/22
     fprintf(fid,'%2d      [%17.8e,%17.8e]   [%17.8e,%17.8e]\n',e,inf(estress(e)),sup(estress(e)),inf(estrain(e)),sup(estrain(e)));
end
%extracting stresses in elements
if (strainflag == 1)
    strain=zeros(ne,1)*infsup(1.,1.);
    j=0;
    for i=ndofnew-ne+1:ndofnew
    j=j+1;
    strain(j) = us(i); % changed u1 to u1s  RLM 9/8/22
    end
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

      fprintf(fid,'%3d       -----              [%15.6e,%15.6e]         ----         ----        -----\n',e,inf(strain(e)),sup(strain(e)));
  end
end %of for loop on elements 

end
return
end

