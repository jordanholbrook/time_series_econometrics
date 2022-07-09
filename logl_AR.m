function L = logl_AR(b)

% The following constructs the loglikelihood function for an AR(1).

global y T AR 

b1 = b(1);                                                                  % Intercept.
b2 = b(2);                                                                  % AR coefficient.
s2 = 
s2(1)  = b(3);                                                                  % Standard deviation.

b1=.5
b2 = .5
s2=1
for t = 2:T
    s2(t) = omega + beta*s2(t-1)+ alpha*(u(t-1))^2
end
L = - 0.5*log(s^2/(1-b2^2)) - 0.5*((y(1) - b1/(1-b2))^2/(s^2/(1-b2^2)));

for t = 2:T
    L = L - 0.5*log(s^2) - 0.5*(((y(t)-b1-b2*y(t-1))^2/s^2));
end

L = L - 0.5*T*log(2*pi);
L = -L;                                                                     % Negative of loglikelihood function (for minimization).

end

