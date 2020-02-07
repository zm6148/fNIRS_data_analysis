function [t, p] = t_test( A, B )
%t-test at each time point

S_A2 = var(A);
S_B2 = var(B);

n_A = length(A);
n_B = length(B);

df = n_A + n_B -2;

S_AB = sqrt(((n_A-1) * S_A2 + (n_B-1) * S_B2)/(n_A + n_B -2));

t = (mean(A) - mean(B))/(S_AB * sqrt(1/n_A + 1/n_B));

p =  betainc(df./(df + t.^2), df/2, 0.5);
end

