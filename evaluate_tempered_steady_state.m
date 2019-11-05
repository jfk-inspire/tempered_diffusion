function y = evaluate_tempered_steady_state(x,alpha,lambda)
% y = evaluate_tempered_steady_state(x,alpha,lambda)
% x = vector of input values
% alpha = fractional order, 1 < alpha <=2
% lambda = tempering parameter, lambda > 0

if (alpha > 2 || alpha <= 1)
    error('alpha must be in interval (1,2]');
end
if (lambda < 0)
    error('lambda must be non-negative')
end
L = min(x);
R = max(x);

if (L >= R)
    error('L must be strictly less than R')
end



diam = R-L;
mid = (R+L)/2;
x1 = (x - L)./diam;


alpham1 = alpha -1;
alpham2 = alpha -2;

if (lambda > 0)

fac1 = mlf_star(alpha,alpham2,lambda.*x1);
fac2 = exp(-lambda.*x1);
u_steady1 = fac1 .* fac2;

fac1 = mlf_star(alpha,alpham1,lambda.*x1);
fac2 = exp(-lambda.*x1);
u_steady2 = fac1 .* fac2;

u_steady = -u_steady2 + u_steady1;

%Normalization constant
%Aconst = (lambda^alpham2)*exp(-lambda*diam)*(diam^alpham1)*mlf(alpha,alpha,(lambda*diam)^alpha);
Aconst = (lambda^alpham2)*exp(-lambda)*mlf(alpha,alpha,lambda^alpha);
y = u_steady./Aconst./diam;

else
    y = alpham1.*x1.^(alpha-2);
    y = y./diam;
end


end



function y = mlf_star(alpha,beta,x)
% Evaluates modified Mittag-Leffler function E^{*}_{a, b} (x)
% See page 89 of Harish's Ph.D. thesis

arg = x.^alpha;
y = x.^beta .* mlf(alpha,beta + 1,arg);

end
