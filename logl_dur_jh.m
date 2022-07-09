function L = logl_dur_jh(b0)

% The following constructs the loglikelihood function for an AR(1).

global x t T N

%b0 = [0 0.5];
XB = [ones(N,1) x]*b0';

  

L(t<40) = log(exp(XB(t<40)))-exp(XB(t<40).*t(t<40));
L(t>=40) = -exp(XB(t>=40)).*40;
L =-sum(L);                                                                % Negative of loglikelihood function (for minimization).

end