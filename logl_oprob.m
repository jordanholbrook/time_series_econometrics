function L = logl_oprob(b0)
global x z 

L = (z==0).*log(normcdf(b0(2)-b0(1).*x)) + (z==1).*log(normcdf(b0(3)-b0(1).*x)-normcdf(b0(2)-b0(1).*x)) + (z==2).*log(1-normcdf(b0(2)-b0(1).*x));

L = -sum(L);
end

