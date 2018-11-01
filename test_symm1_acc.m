function power_symm = test_symm1_acc(M,N,setting,n,p)
N; % random permutations
M; % simulations
T_symm = []; % Symmetry Test Statistics
pT_symm = []; % Symmetry Permutation Test Statistics
powers = [];
h = (5 - (-5))/(p-1); mu = [-5:h:5];
for m = 1:M
  m;
    switch setting
      %---------------------------------Symmetry distribution---------------------------------%
      case 1 % Multivariate Normal distribution (dependent structure 1)
        sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
        X = mvnrnd(mu,sigma,n);
      case 2 % Multivariate Normal distribution (dependent structure 2)
        sigma = diag(ones(1,p));
        c_max = 0.9; c_min = 0.1;
        delta = (c_max - c_min)/(nchoosek(p,2)-1);
        c = c_min:delta:c_max;
        idx = tril(true(size(sigma)), -1); % Lower triangular half
        sigma(idx) = c;
        sigma = tril(sigma,-1)' + sigma;
        X = mvnrnd(mu,sigma,n);
      case 3 % Mixed Multivariate Normal distribution (dependent structure 1&2)
        sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1; X1 = mvnrnd(mu,sigma,0.5*n);
        sigma = diag(ones(1,p));
        c_max = 0.9; c_min = 0.1;
        delta = (c_max - c_min)/(nchoosek(p,2)-1);
        c = c_min:delta:c_max;
        idx = tril(true(size(sigma)), -1); % Lower triangular half
        sigma(idx) = c;
        sigma = tril(sigma,-1)' + sigma;
        X2 = mvnrnd(mu,sigma,0.5*n);
        X = [X1;X2];
      case 4 % Mixed Multivariate Normal distribution (dependent structure 1)
      sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1; X1 = mvnrnd(mu,sigma,0.9*n);
      sigma = diag(ones(1,p));
      c_max = 0.9; c_min = 0.1;
      delta = (c_max - c_min)/(nchoosek(p,2)-1);
      c = c_min:delta:c_max;
      idx = tril(true(size(sigma)), -1); % Lower triangular half
      sigma(idx) = c;
      sigma = tril(sigma,-1)' + sigma;
      X2 = mvnrnd(mu,sigma,0.1*n);
      X = [X1;X2];
      case 5 % Mixed Multivariate t distribution
        v = 3; sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
        X = ((sigma.^0.5)*(trnd(v,n,p))')';
      %-------------------------------alternative distribution--------------------------------%
    case 6 % chi^2
        v = 3; sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
        temp = chi2rnd(v,n,p)- v + repmat(mu,[n,1]);
        X = ((sigma.^0.5)*(temp)')';
      case 7 % chi^2
        v = 5; sigma = eye(p,p);
        X = chi2rnd(v,n,p)- v + repmat(mu,[n,1]);
        % X = ((sigma.^0.5)*(temp)')';
      case 8 % chi^2
        v = 6; sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
        temp = chi2rnd(v,n,p)- v + repmat(mu,[n,1]);
        X = ((sigma.^0.5)*(temp)')';
      case 9 % Mixed Multivariate Normal distribution and ��^2
        sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1; X1 = mvnrnd(mu,sigma,0.9*n);
        v = 4; sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
        temp = chi2rnd(v,0.1*n,p)- v + repmat(mu,[0.1*n,1]);
        X2 = ((sigma.^0.5)*(temp)')';
        X = [X1;X2];
      case 10 % Gamma distribution
        k = 2; % shape parameter [��(1;2;3;5;10)]
        theta = 2; % scale parameter
        sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
        temp = gamrnd(k,theta,n,p)- k*theta + repmat(mu,[n,1]);
        X = ((sigma.^0.5)*(temp)')';
      case 11 % Gamma distribution
        k = 3; % shape parameter [��(1;2;3;5;10)]
        theta = 2; % scale parameter
        sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
        temp = gamrnd(k,theta,n,p)- k*theta + repmat(mu,[n,1]);
        X = ((sigma.^0.5)*(temp)')';
      case 12 % Gamma distribution
        k = 5; % shape parameter [��(1;2;3;5;10)]
        theta = 2; % scale parameter
        sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
        temp = gamrnd(k,theta,n,p)- k*theta + repmat(mu,[n,1]);
        X = ((sigma.^0.5)*(temp)')';
      case 13 % Mixed Multivariate t distribution
       v = 5; sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
       X = ((sigma.^0.5)*(trnd(v,n,p))')';
      case 14 % Gamma distribution
       k = 1; % shape parameter [��(1;2;3;5;10)]
       theta = 2; % scale parameter
       sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
       temp = gamrnd(k,theta,n,p)- k*theta + repmat(mu,[n,1]);
       X = ((sigma.^0.5)*(temp)')';
     case 15 % chi^2
       v = 3; sigma = 0.5*ones(p,p); sigma(1:p+1:end) = 1;
       temp = chi2rnd(v,n,p)- v + repmat(mu,[n,1]);
       X = ((sigma.^0.5)*(temp)')';
      %-------------------------------other distribution--------------------------------%
      case 16 % multivariate standard normal distribution
        mu_16 = zeros(1,p);
        sigma = eye(p,p);
        X = mvnrnd(mu_16,sigma,n);

      end

    Y = X - repmat(mean(X),[n,1]);
    T_symm = symm(Y)';

    % Generate U with P(U=1)=0.5��P(U=-1)=0.5;
    U_temp = binornd(1,0.5,n*N,1);
    U = (2*U_temp - 1) * ones(1,p);
    U = permute(reshape(U',[p,n,N]), [2,1,3]);
    %----------------special bootstrap-----------------%
    Y_boot_temp = repmat(Y,[1,1,N]).*U;
    Y_boot = Y_boot_temp - repmat(mean(Y_boot_temp),[n,1,1]);
    pT_symm = symm_boot(Y_boot);

    % calculate power
    pvalue = mean(pT_symm >= repmat(T_symm,[1,N]),2);
    powers(m,:) = (pvalue < 0.05);
end
power_symm = mean(powers);
end

% compute the test statistic corresponding to the Symmetry Test Statistic
% input: Y      nxp     independent sample matrix
%        a      (0,2]   characteristic exponent
% output:T      observed test statistic
function T = symm(Y)
n = size(Y,1); % sample size
Y_subtr_norm = [];
Y_add_norm   = [];
Y_subtr_norm = sum((kron(Y,ones(n,1)) - kron(ones(n,1),Y)).^2,2);
Y_add_norm   = sum((kron(Y,ones(n,1)) + kron(ones(n,1),Y)).^2,2);
T = [];
c = 1;
for a = 0.5:0.5:2.0
T(c) = (1/n^2) * sum( exp(-Y_subtr_norm.^(a/2)) - exp(-Y_add_norm.^(a/2))  );
c = c + 1;
end
end


function T = symm_boot(Y_boot)
n = size(Y_boot,1); % sample size
p = size(Y_boot,2); % sample dimension
N = size(Y_boot,3); % bootstrap size
Y_boot_subtr_norm = [];
Y_boot_add_norm   = [];
T = [];
% using matrix transform to perform 3D kron
Y_boot_subtr_norm = sum((permute(reshape(repmat(permute(Y_boot,[2,1,3]),[n,1,1]),[p,n*n,N]),[2,1,3]) - repmat(Y_boot,[n,1,1])).^2,2);
Y_boot_add_norm   = sum((permute(reshape(repmat(permute(Y_boot,[2,1,3]),[n,1,1]),[p,n*n,N]),[2,1,3]) + repmat(Y_boot,[n,1,1])).^2,2);
c = 1;
for a = 0.5:0.5:2.0
T(c,:) = (1/n^2) * sum( exp(-Y_boot_subtr_norm.^(a/2)) - exp(-Y_boot_add_norm.^(a/2))  );
c = c + 1;
end
end
