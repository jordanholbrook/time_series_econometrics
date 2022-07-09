function L = logl_order_prob(b0)

% The following constructs the loglikelihood function for an ordered probit model.

global x z 


b1=b0(1);

XB = b1*x;

L(z == 0) = log(normcdf(b0(2)- XB(z==0)));
L(z == 1) = log(normcdf(b0(3)- XB(z==1)) -normcdf(b0(2)- XB(z==1)));
L(z == 2) = log(1 - (normcdf(b0(3)-XB(z==2))));
L = -sum(L);
end