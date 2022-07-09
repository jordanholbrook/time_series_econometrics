function L = logl_MA2( b )

% The following constructs the loglikelihood function for an MA(2).

global x2 T 

omega = zeros(T,T);                                                         % Placeholder for Variance-Covariance Matrix.
b0 = b(1);                                                                  % Mean.
b1 = b(2);
b2 = b(3); 
s  = b(4);                                                                  % Standard deviation.

omega = omega + eye(T).*(s^2).*(1+b1^2+b2^2);                                    % Fill in the variances.
        
omega(2,1) = (s^2)*(b1+b1*b2);                                                      % Fill in the first order covariances.
omega(T-1,T) = omega(2,1);
omega(3,1) = (s^2)*b2;
omega(T-2,T) = omega(3,1);

for i = 2:T-1
    omega(i-1,i) = (s^2)*(b1+b1*b2);
    omega(i+1,i) = omega(i-1,i);
end    

for i = 3:T-2
    omega(i-2,i) = (s^2)*b2;
    omega(i+2,i) = omega(i-2,i);
end
    
L = - 0.5*log(abs(det(omega)))- 0.5*(x2-b0)'*inv(omega)*(x2-b0);             
L = L - 0.5*T*log(2*pi);
L = -L;                                                                     % Negative of loglikelihood function (for minimization).

end

