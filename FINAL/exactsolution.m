
function fn = exactsolution(x)

if(x>=0.6 && x <=0.8)
    fn = 1.0;
elseif(x>=0.2 && x<=0.4)
    fn = sin(5.0*pi*(x-0.2));
else 
    fn = 0.0; 
end
