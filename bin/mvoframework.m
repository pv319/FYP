function [max_sharpe_weights, max_sharpe_returns, max_sharpe_risk] = mvoframework(mean_returns, cov_matrix, shortselltoggle, rf)


% %excel file import
% filename = 'fulldata.xlsx';
% cdata = readmatrix(filename);
% disp(size(cdata));
% 
% cdata1 = cdata(:,:);
% 
% [n, m] = size(cdata1);
% returns = zeros(n-1, m);
% for j = 1:m
%     for i = 2:n
%         returns(i-1, j) = (cdata1(i, j) - cdata1(i-1, j)) / cdata1(i-1, j);
%     end
% end

%%finding expected ret, sigma and the size
% mean_returns = mean(returns)';
% cov_matrix = cov(returns);
num_assets = size(mean_returns, 1);

% rf rate for sharpe calculation
rf = 0.00;

% we aim to find the lowest risk possible for each interval between
% returns. Accuracy increases with the more number of portfolios.
target_returns = linspace(min(mean_returns), max(mean_returns), 100);

% Array setup
portfolio_risks = zeros(length(target_returns), 1);
portfolio_returns = zeros(length(target_returns), 1);
portfolio_weights = zeros(length(target_returns), num_assets);

% Solve the optimization problem for each linspaced return
options = optimoptions('quadprog', 'Display', 'off');
for i = 1:length(target_returns)
    target_return = target_returns(i);
    
    % Define the optimization problem
    f = zeros(num_assets, 1); %no constant
    Aeq = ones(1, num_assets); %sum(weights) ==1, we can change this to -1 to 1 to allow for shortselling
    beq = 1;
    A = -mean_returns'; % mean return >= target return
    b = -target_return;
    if (shortselltoggle == 0    )
        lb = zeros(num_assets, 1);
        ub = ones(num_assets, 1);
    else
        lb = -1 *ones(num_assets, 1);
        ub = ones(num_assets, 1);
    end

    % Solve the quadratic programming problem
    [x, ~, solutionindicator] = quadprog(cov_matrix, f, A, b, Aeq, beq, lb, ub, [], options);
    
    if solutionindicator == 1 % Solution found
        portfolio_weights(i, :) = x';
        portfolio_returns(i) = mean_returns' * x;
        portfolio_risks(i) = sqrt(x' * cov_matrix * x);
    else
        portfolio_weights(i, :) = NaN;
        portfolio_returns(i) = NaN;
        portfolio_risks(i) = NaN;
    end
end

% Remove NaN values
valid_indices = ~isnan(portfolio_returns);
portfolio_risks = portfolio_risks(valid_indices);
portfolio_returns = portfolio_returns(valid_indices);
portfolio_weights = portfolio_weights(valid_indices, :);

%find the sharpe ratios for each generated portfolio
sharpe_ratios = (portfolio_returns - rf) ./ portfolio_risks;

%locate portfolio with max sharpe ratio
[~, sharpehigh] = max(sharpe_ratios);
%return parameters
max_sharpe_weights = portfolio_weights(sharpehigh, :);
max_sharpe_return = portfolio_returns(sharpehigh);
max_sharpe_risk = portfolio_risks(sharpehigh);



%Efficient Frontier plot with risk free rate line (CML)
figure;
hold on;
plot(portfolio_risks, portfolio_returns, 'b-', 'LineWidth', 2);

plot(max_sharpe_risk, max_sharpe_return, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Risk (Standard Deviation)');
ylabel('Return');
title('Efficient Frontier and Maximum Sharpe Ratio Portfolio');
legend('Efficient Frontier', 'Max Sharpe Ratio Portfolio', 'Location', 'Best');
grid on;

%
disp('weights for the high sharpe ratio portfolio');
disp(max_sharpe_weights);

disp(['Maximum Sharpe Ratio: ', num2str(sharpe_ratios(sharpehigh))]);
disp(['Expected Portfolio Return: ', num2str(max_sharpe_return)]);
disp(['Portfolio Risk (Standard Deviation): ', num2str(max_sharpe_risk)]);
end