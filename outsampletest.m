function outsampletest(pwgt1, outsamples)

muoutsample = mean(outsamples);
disp('avg returns out sample');
disp(muoutsample');
Sigmaoutsample = cov(outsamples);

%mkmvo for outsample data
% Create a portfolio object
p = Portfolio;

% Set asset mean returns and covariance matrix
p = p.setAssetMoments(muoutsample, Sigmaoutsample);

% Set portfolio constraints
p = p.setDefaultConstraints;

% Perform mean-variance portfolio optimization
pwgt2 = estimateFrontier(p);

%insample weights
disp(pwgt2);
%outsample weights
portfolio_stdev = zeros(length(pwgt1),1);
portfolio_stdev2 = zeros(length(pwgt2),1);
for i = 1:length(pwgt1)
    portfolio_stdev(i) = sqrt(pwgt1(:,i)' * Sigmaoutsample * pwgt1(:,i));
    portfolio_stdev2(i) = sqrt(pwgt2(:,i)' * Sigmaoutsample * pwgt2(:,i));
end
portfolio_returns1 = pwgt1' * muoutsample';
portfolio_returns2 = pwgt2'*muoutsample';

maxmin = 0;
maxminpf = 0;
temp = 0;
for j = 1:length(pwgt2)
    temp = portfolio_returns2(j) / portfolio_stdev2(j);
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

plot(portfolio_stdev2, portfolio_returns2, 'x', 'MarkerSize', 10, 'LineWidth', 2);

% Label the points
for i = 1:length(pwgt1)
    if (i == maxminpf)
        
        text(portfolio_stdev(maxminpf), portfolio_returns1(maxminpf), sprintf('Portfolio %d, most optimal portfolio', i),'Color','r','HorizontalAlignment','left','VerticalAlignment','bottom');
        i = i+1;
    end
    %text(portfolio_stdev(i), portfolio_returns1(i), sprintf('Portfolio %d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

%text(portfolio_stdev(maxminpf), portfolio_returns1(maxminpf), 'Best Portfolio', 'FontSize',12,'Color','r','HorizontalAlignment','right','VerticalAlignment','bottom');

xlabel('Standard Deviation');
ylabel('Expected Return');
title('Efficient Frontier with Optimal Portfolio Weights taken from insample tested on returns of outsample');
hold off;

optimalweight = pwgt1(:,maxminpf);
disp('the optimal weightage for the 10 assets is');
disp(optimalweight);

%Calculate the Correlation between the two datasets.
correlation_matrix = corrcoef([portfolio_stdev, portfolio_returns1], [portfolio_stdev2, portfolio_returns2]);
correlation_coefficient = correlation_matrix(1, 2);
fprintf('Correlation coefficient between the two groups: %.4f\n', correlation_coefficient);

%Calculate the euclidean distance between the two datasets
edist = 0;
etemp = 0;
for i = 1:length(pwgt1)
    etemp = sqrt(((portfolio_stdev2(i)-portfolio_stdev(i))^2) - ((portfolio_returns2(i)-portfolio_returns1(i))^2));
    edist = edist + etemp;
end
ecliddistavg = edist/length(pwgt1);

printout = sprintf('euclidean distance between pwgt1 and pwgt2 is %d', ecliddistavg);
disp(printout);

initial_value = 1;
portfolio_value = initial_value * cumprod(1 + (outsamples * pwgt1(5)));

% figure;
% plot(portfolio_value);
% xlabel('Time');
% ylabel('Portfolio Value ($)');
% title('Portfolio Value Over Time (Out-of-Sample)');
% grid on;




end