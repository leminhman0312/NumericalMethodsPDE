%exact solution 

function fn = uexact(t,x,y)
    fn = exp(-t) sin(pi*x) sin(pi*y)
end
