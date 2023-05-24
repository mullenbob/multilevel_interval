function [val] = dedtobin(decimal,m)
i=0;
while(decimal > 0)
   i=i+1;
   rem(i) = decimal -2*fix(decimal / 2);
   decimal = fix(decimal/2);
end
n = i;
m = max(n,m);
val = zeros(m,1);
for i=1:m-n
   val(i,1)=0;
end
for i=m-n+1:1:m
  val(i,1)= rem(m-i+1);
end
val = val';
return
