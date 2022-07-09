function L = logl_MA( b )

% The following constructs the loglikelihood function for an MA(1).

global x T MA

omega = zeros(T,T);                                                         % Placeholder for Variance-Covariance Matrix.
b0 = b(1);                                                                  % Mean.
b1 = b(2);                                                                  % MA coefficient.
s  = b(3);                                                                  % Standard deviation.

omega = omega + eye(T).*(s^2).*(1+b1^2);                                    % Fill in the variances.
        
omega(2,1) = (s^2)*b1;                                                      % Fill in the first order covariances.
omega(T-1,T) = omega(2,1);

for i = 2:T-1
    omega(i-1,i) = (s^2)*b1;
    omega(i+1,i) = omega(i-1,i);
end    
    
L = - 0.5*log(abs(det(omega)))- 0.5*(x-b0)'*inv(omega)*(x-b0);              % Loglikelihood function.
L = L - 0.5*T*log(2*pi);
L = -L;                                                                     % Negative of loglikelihood function (for minimization).

end

