function [sigma_shrinkage, delta] = ledoitWolf(returns, manualdelta)
    % Implementation of the Ledoit-Wolf shrinkage estimator for covariance matrix
    
    
    [n, p] = size(returns);
    
    
    empirical_cov = cov(returns);
    
    % Mean of the empirical covariance matrix diagonal
    mean_variance = mean(diag(empirical_cov));
    
    % Structured estimator (diagonal matrix with the mean variance on the diagonal)
    structured_cov = mean_variance * eye(p);
    
    %     # sample average correlation
    % var = np.diag(sample_cov).reshape(-1, 1)
    % sqrt_var = var ** 0.5
    % unit_cor_var = sqrt_var * sqrt_var.transpose()
    % average_cor = ((sample_cov / unit_cor_var).sum() - n) / n / (n - 1)
    % prior = average_cor * unit_cor_var
    % np.fill_diagonal(prior, var)
    % 
    % # pi-hat
    % y = returns ** 2
    % phi_mat = (y.transpose() @ y) / t - sample_cov ** 2
    % phi = phi_mat.sum()
    % 
    % # rho-hat
    % theta_mat = ((returns ** 3).transpose() @ returns) / t - var * sample_cov
    % np.fill_diagonal(theta_mat, 0)
    % rho = (
    %     np.diag(phi_mat).sum()
    %     + average_cor * (1 / sqrt_var @ sqrt_var.transpose() * theta_mat).sum()
    % )

    % phi = 0;
    % for i = 1:n
    %     diff = returns(i, :)' * returns(i, :) - empirical_cov;
    %     phi = phi + norm(diff, 'fro')^2;
    % end
    % phi = phi / n;
    % 
    % 
    % rho = 0;
    % for i = 1:n
    %     diff = returns(i, :)' * returns(i, :) - structured_cov;
    %     rho = rho + norm(diff, 'fro')^2;
    % end
    % rho = rho / n;
    % 
    % 
    % gamma = norm(empirical_cov - structured_cov, 'fro')^2;
    % 
    % 
    % kappa = (phi - rho) / gamma;
    % delta = max(0, min(1, kappa / n));
    % if (delta == 0)
    %     delta = manualdelta;
    %     disp('delta was calculated as 0, therefore manualdelta has been usen for shrinkage factor');
    % end
    
    % Shrinkage covariance matrix
    sigma_shrinkage = delta * structured_cov + (1 - delta) * empirical_cov;
end