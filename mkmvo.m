function [pwgt1,outSampleReturns] = mkmvo(cdata, k)

% period = 1;
% range1 = 1 + (24*period);
% range2 = (period+1)*12;
cdata1 = cdata(:,:);

% Calculate returns
%returns = diff(cdata1) ./ cdata1(1:end-1,:);

[n, m] = size(cdata1);
returns = zeros(n-1, m);
for j = 1:m
    for i = 2:n
        returns(i-1, j) = (cdata1(i, j) - cdata1(i-1, j)) / cdata1(i-1, j);
    end
end

% Divide data into in-sample and out-of-sample sets
inSampleSize = round(0.5 * size(returns, 1)); % 70% in-sample
inSampleSizeData = round(0.5 * size(cdata1, 1));
inSampleReturns = returns(1:inSampleSize, :);
outSampleReturns = returns(inSampleSize+1:end, :);
outSampleData = cdata1(inSampleSizeData+1:end, :);

% Define expected return and covariance matrix
mu = mean(inSampleReturns);
Sigma = cov(inSampleReturns);
disp(mu');


% Create a portfolio object
p = Portfolio;
numportfolios = 50;

% Set asset mean returns and covariance matrix
p = p.setAssetMoments(mu, Sigma);

% Set portfolio constraints
p = p.setDefaultConstraints;
% Set portfolio constraints to allow short selling
% Allow weights to be between -1 and 1
%p = p.setBounds(-1, 1);
% p = p.setNumPorts(100);

% Perform mean-variance portfolio optimization
pwgt1 = estimateFrontier(p, numportfolios);
sharpe = estimateMaxSharpeRatio(p);
disp('sharpe');
disp(sharpe);
%a = pwgt1(:,1);

% Display optimal portfolio weights
disp(pwgt1);
portfolio_stdev = zeros(length(pwgt1),1);
for i = 1:length(pwgt1)
    portfolio_stdev(i) = sqrt(pwgt1(:,i)' * Sigma * pwgt1(:,i));
end
portfolio_returns1 = pwgt1' * mu';
temp = 0;
maxminpf = 1;
maxmin = portfolio_returns1(1)/((k)*portfolio_stdev(1));
for j = 1:length(pwgt1)
    temp = portfolio_returns1(j)/((k)*portfolio_stdev(j));
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
for i = 1:length(pwgt1)
    if (i == maxminpf)
        text(portfolio_stdev(maxminpf), portfolio_returns1(maxminpf), sprintf('Portfolio %d, most optimal portfolio', i),'Color','r','HorizontalAlignment','left','VerticalAlignment','bottom');
        
    end
    text(portfolio_stdev(i), portfolio_returns1(i), sprintf('Portfolio %d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    i = i+1;
end

%text(portfolio_stdev(maxminpf), portfolio_returns1(maxminpf), 'Best Portfolio', 'FontSize',12,'Color','r','HorizontalAlignment','right','VerticalAlignment','bottom');

xlabel('Standard Deviation');
ylabel('Expected Return');
title('Efficient Frontier with Optimal Portfolio Weights (markowitz insample)');
hold off;

optimalweight = pwgt1(:,maxminpf);
disp('the optimal weightage for the 10 assets is');
disp(optimalweight);


% initial_value = 1;
% portfolio_value = 0;
% %portfolio_value = initial_value * cumprod(1 + (outSampleData * pwgt1(maxminpf,:)'));
% for i = 1:length(outSampleData)
% 
% 
% figure;
% plot(portfolio_value);
% xlabel('Time');
% ylabel('Portfolio Value ($)');
% title('Portfolio Value Over Time (Out-of-Sample)');
% grid on;
% 

end