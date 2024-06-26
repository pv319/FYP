function [dateMatrix, stockMatrix] = datafetch(tickers, apiKey)
    % fetchStockData Fetches historical daily closing stock data for a set of tickers from Alpha Vantage
    % for the date range from January 1, 2013 to January 1, 2023, and stores the data in an n x m matrix.
    %
    % Parameters:
    %   tickers (cell array): A cell array of stock tickers.
    %   apiKey (char): Your Alpha Vantage API key.
    %
    % Returns:
    %   dateMatrix (datetime array): Array of dates.
    %   stockMatrix (double array): n x m matrix with closing prices where n is the number of dates
    %                               and m is the number of companies.

    baseUrl = 'https://www.alphavantage.co/query';
    stockData = struct();

    startDate = datetime('2013-01-01');
    endDate = datetime('2023-01-01');
    numTickers = length(tickers);

    % Initialize the date matrix
    allDates = (startDate:endDate)';
    numDates = length(allDates);

    % Initialize the stock matrix with NaNs
    stockMatrix = NaN(numDates, numTickers);

    % Iterate over each ticker
    for i = 1:numTickers
        ticker = tickers{i};
        params = {
            'function', 'TIME_SERIES_DAILY', ...
            'symbol', ticker, ...
            'outputsize', 'full', ...
            'apikey', apiKey
        };

        response = webread(baseUrl, params{:});
        if isfield(response, 'TimeSeries_Daily_')
            data = response.TimeSeries_Daily_;

            % Convert the nested structure to a table
            data = struct2table(data);
            data = table2array(data);

            % Transpose the array
            data = data';

            % Convert the transposed array back to a table
            data = array2table(data);
            data.Properties.VariableNames = strrep(strrep(data.Properties.VariableNames, 'x1_', ''), '_', '');

            % Check if the 'close' column exists
            if ismember('close', data.Properties.VariableNames)
                % Extract only the 'close' data and filter by date range
                dates = datetime(data.Properties.RowNames, 'InputFormat', 'yyyy-MM-dd');
                closeData = data(:, {'close'});
                closeData.Date = dates;

                % Filter the data between the specified date range
                closeData = closeData(closeData.Date >= startDate & closeData.Date <= endDate, :);

                % Align the close prices with the allDates array
                [~, idx] = ismember(closeData.Date, allDates);
                stockMatrix(idx, i) = str2double(closeData.close);
            else
                warning('No closing price data found for ticker: %s', ticker);
            end
        else
            warning('No data found for ticker: %s', ticker);
            % Print the raw response for debugging purposes
            disp(['Raw response for ticker ', ticker, ':']);
            disp(response);
        end
    end

    % Output the date matrix
    dateMatrix = allDates;
end