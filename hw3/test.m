clear all 
clc

x = 1:10;
n = length(x);
avg = mymean(x,n);

myprint(2);


k = myprint(2);

disp(k);

function a = mymean(v,n)
% MYMEAN Example of a local function.

    a = sum(v)/n;
end




function h = myprint(x)
  
  x = x+10;
  h = x;
  
end