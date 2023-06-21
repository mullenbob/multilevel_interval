function secant= secant(kode,strain,E,Eplastic)
%Compute the interval secant modulus for nonlinear materials
% kode defines the model  1=cubic  2= %Ramberg Osgood model   n= Eplastic S1=fYR
% kode = 3 bi linear  

if (kode == 1) 
   %%cubic  stress strain model  E * strain + Eplastic*(strain^3)
   %upper bound of strain 
   supstrain=sup(strain(e));
   % fix for very small strain
   if (abs(supstrain) < 1.E-50)
       maxstress=0.;
       secant1=sup(E);
   else
  %check on weither we should subtract infstrain here dependency issue  RLM 5/26/2023
    maxstress=sup(E)*supstrain-Eplastic*supstrain^3;
   
   secant1=maxstress/supstrain;
   end
   infstrain=inf(strain(e));
   if (abs(infstrain) < 1.e-50) 
       minstress=0.;
       secant2=inf(E);
   else
       minstress=inf(E)*infstrain-Eplastic*infstrain^3;
   
       secant2=minstress/infstrain;
   end   
    
end
if (kode == 2) 
   %Ramberg Osgood model   n= Eplastic S1=fY
   supstrain=sup(strain(e));
   if (abs(supstrain) < 1.E-50)
       maxstress=0.;
       secant1=E;
   else
   [maxstress  tan]=ramberg(supstrain,E,fY,Eplastic);    
   
   secant1=maxstress/supstrain;
   end
   infstrain=inf(strain(e));
   if (abs(infstrain) < 1.e-50) 
       minstress=0.;
       secant2=E;
   else
        [minstress,  tan]=ramberg(infstrain,E,fY,Eplastic);
        secant2=minstress/infstrain;
   end  
    
end

        
secant = infsup(min(secant1,secant2),max(secant1,secant2));
return
end
function [stress, tan]= ramberg(strain,E,S1,N)
tol=10E-6;
%program to calculate the Ramberg-Osgood one dimensionial consititutive model
%%A PROGRAM BY Dr. RAMA RAO MALLELA , Prof. RAFI MUHANNA AND Prof.  ROBERT MULLEN
% strain = stress/E + K*(stress/E)^n   using newton method to solve
%use elastic value as initial guess
%convert to non-dimensional version witn S1 and m1=.7 from paper 
%NACA  TN-902  1943
strain1=strain*E/S1;
if (strain > 0) 
stress1=strain1;
delta=1;
count=1;
while (abs(delta) > tol)
% for i=1:3
count=count+1;
if (count>1000) 
    break;
end
     delta=strain1-stress1-3/7*(stress1)^N;
     tan=(1.+3./7*N*stress1^(N-1));
     
stress1=stress1+delta/tan;
 end
stress=stress1*S1;
tan=E/tan;
else
    strain1=-strain1;
    stress1=strain1;
delta=1;
while (abs(delta) > tol)

     delta=strain1-stress1-3/7*(stress1)^N;
     tan=(1.+3./7*N*stress1^(N-1));
     
stress1=stress1+delta/tan;
 end
stress=-stress1*S1;
tan=E/tan;
end
return
end