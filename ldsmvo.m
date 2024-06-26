function [pwgt, out_sample_returns] = ldsmvo(cdata, manualdelta)
cdata1 = cdata(:,:);

% Calculate returns
returns = diff(cdata1) ./ cdata1(1:end-1,:);

%manula delta mmust vary from 0-1
%manualdelta = 0;

% Split into in-sample and out-of-sample
inSampleSize = round(0.5 * size(returns, 1)); % 70% in-sample
in_sample_returns = returns(1:inSampleSize, :);
out_sample_returns = returns(inSampleSize+1:end, :);

% Use Ledoit-Wolf shrinkage to estimate the covariance matrix
[Sigma_LW, delta] = ledoitWolf(in_sample_returns, manualdelta);

% % Display the shrinkage intensity parameter
disp(['Shrinkage Intensity Parameter (delta): ', num2str(delta)]);


% Create a portfolio object
p = Portfolio;

mu_LW = mean(in_sample_returns);

% Set asset mean returns and covariance matrix
p = p.setAssetMoments(mu_LW, Sigma_LW);

% Set portfolio constraints
p = p.setDefaultConstraints;

% Perform mean-variance portfolio optimization
pwgt = estimateFrontier(p);
a = pwgt(:,1);

% Display optimal portfolio weights
disp(pwgt);
portfolio_stdev = zeros(length(pwgt),1);
for i = 1:length(pwgt)
    portfolio_stdev(i) = sqrt(pwgt(:,i)' * Sigma_LW * pwgt(:,i));
end
portfolio_returns1 = pwgt' * mu_LW';
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