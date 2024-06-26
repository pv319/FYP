 
function [pwgt, out_sample_returns] = ssmvo(cdata, evc)


cdata1 = cdata(1:2517,:);


returns = diff(cdata1) ./ cdata1(1:end-1, :);

% Split into in-sample and out-of-sample
inSampleSize = round(0.5 * size(returns, 1)); % 70% in-sample
in_sample_returns = returns(1:inSampleSize, :);
out_sample_returns = returns(inSampleSize+1:end, :);
% in_sample_returns = returns(1:15, :); % First 15 observations for in-sample
% out_sample_returns = returns(16:end, :); % Remaining for out-of-sample


in_sample_mean = mean(in_sample_returns)';
in_sample_cov = cov(in_sample_returns);

% doing pca
[coeff, latent, explained] = pcacov(in_sample_cov);

%
explained_variance = cumsum(explained);
num_components = find(explained_variance >= evc, 1); %we can change this, varying this allows us to experiment and decide what we can change here, how it affects computation time
disp(num_components);

% dimensionality reduction
coeff_reduced = coeff(:, 1:num_components);
latent_reduced = latent(1:num_components); %check this

% mean is smaller subsace, variance too
%have to project via transpose of the new coeffs
mean_reduced = coeff_reduced' * in_sample_mean;
cov_reduced = diag(latent_reduced);
%cov_reduced = coeff_reduced' * latent_reduced * coeff_reduced;

% bring back to original space
in_sample_mean_reconstructed = coeff_reduced * mean_reduced;
%to create pxn matrix, multiple x (diagonal)  x(^t)
in_sample_cov_reconstructed = coeff_reduced * cov_reduced * coeff_reduced';

% % Create a portfolio object
p = Portfolio;


% Set asset mean returns and covariance matrix
p = p.setAssetMoments(in_sample_mean_reconstructed, in_sample_cov_reconstructed);

% Set portfolio constraints
p = p.setDefaultConstraints;

% Perform mean-variance portfolio optimization
pwgt = estimateFrontier(p);
a = pwgt(:,1);

% Display optimal portfolio weights
disp(pwgt);
portfolio_stdev = zeros(length(pwgt),1);
for i = 1:length(pwgt)
    portfolio_stdev(i) = sqrt(pwgt(:,i)' * in_sample_cov_reconstructed * pwgt(:,i));
end
portfolio_returns1 = pwgt' * in_sample_mean_reconstructed;
maxmin = 0;
maxminpf = 0;
temp = 0;
for j = 1:length(pwgt)
    temp = portfolio_returns1(j) / portfolio_stdev(j);
    disp('current sharpe ratio')
    disp(j);
    disp(temp);
    if (maxmin<temp)
        maxmin = temp;
        maxminpf = j;
    end
end
printout = sprintf('Portfolio %d is the portfolio offering the most return considering risk', maxminpf);
disp(printout);
EFpoints = [portfolio_stdev, portfolio_returns1];
disp(EFpoints);
%coordinates_array = table2array(EFpoints);

% Plot efficient frontier

figure;
hold;
plotFrontier(p);

% Calculate standard deviation for each portfolio

%add risk bias



% Plot optimal portfolio weights on the efficient frontier
plot(portfolio_stdev, portfolio_returns1 , 'o' );

% Label the points
for i = 1:length(pwgt)
    if (i == maxminpf)
        
        text(portfolio_stdev(maxminpf), portfolio_returns1(maxminpf), sprintf('Portfolio %d, most optimal portfolio', i),'Color','r','HorizontalAlignment','left','VerticalAlignment','bottom');
        i = i+1;
    end
    text(portfolio_stdev(i), portfolio_returns1(i), sprintf('Portfolio %d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

%text(portfolio_stdev(maxminpf), portfolio_returns1(maxminpf), 'Best Portfolio', 'FontSize',12,'Color','r','HorizontalAlignment','right','VerticalAlignment','bottom');

xlabel('Standard Deviation');
ylabel('Expected Return');
title('Efficient Frontier with Optimal Portfolio Weights');
hold off;

optimalweight = pwgt(:,maxminpf);
disp('the optimal weightage for the 10 assets is');
disp(optimalweight);

end