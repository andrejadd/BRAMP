

alpha_sigma = 1
beta_sigma = 0.1

i_node = 4;
DAG = [1,1,1,1;0,1,0,1;0,1,0,0];
mu_vec = [1,2,3,4;11,12,13,14;21,22,23,24;31,32,33,34;41,42,43,44];
parents = find(DAG(:,i_node))

mu_vec = mu_vec([1;parents+1],1)

%mu_vec = 0.2
lambda_snr = 0.5

Cov_mat = [1,2,3;1.1,1.8,1.3; 21,27.8,23]

data = [1.1,1.2,1.3,1.4,1.5;
             2.1,2.2,2.3,2.4,2.5;
             3.1,3.2,3.3,3.4,3.5;
             4.1,4.2,4.3,4.4,4.5;  ]

[n_plus, n_obs] = size(data)

X = [ones(1,n_obs);data(parents,:)]

y = data(2,:)'


Sigma_inv = inv(lambda_snr*Cov_mat) + X*X'                        

mue_apost = inv(Sigma_inv)*(inv(lambda_snr*Cov_mat)*mu_vec+X*y)
           
W_i = mvnrnd(mue_apost,var_all*inv(Sigma_inv))
% W_i = W_i'
