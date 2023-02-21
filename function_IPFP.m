function [IPFP] = function_IPFP(val_link, A, TM_estimation)
% The IPFP algorithm is able to adjust the data of the matrix to fit a given constraint through iterative iterations.
% A matrix is represented as a routing matrix, the rows indicate the number of the link, the columns indicate the number of the od flow and whether the od flow passes through the current link
[r_A, c_A] = size(A);
lb = zeros(c_A, 1);

n_link = r_A;                                             % link number
n_flow = c_A;                                             % flow number
n_loop = 70;                                              % loop times
TM = TM_estimation;                                       % Current predicted traffic size matrix

% Set the error bounds used by the IPFP algorithm
y_erro = zeros(n_link, 1);
x_erro = zeros(n_flow, 1);
y_erro(:, 1) = 0.001;                                     % link byte counts error bound
x_erro(:, 1) = 0.001;                                     % OD flow traffic error bound


% The rows represent the total number of od streams and the columns represent the total number of hours
[r_est, c_est] = size(TM);
  
h = waitbar(0, 'please waiting ...');

% Shaping of the entire od stream at each moment, using the IPFP algorithm
for tm = 1:c_est
    waitbar(tm / c_est, h, sprintf('handling traffic matrix at the time %d, please waiting ...', tm));
    
    % Get the per-od stream size and per-link load at the moment of tm
    x_vl2 = TM(:, tm);
    y_cur_vl = val_link(:, tm);

    x_isinf = isinf(x_vl2);              % Determine if the elements in x_vl2 are infinite
    x_inf = find(x_isinf == 1);
    x_vl2(x_inf) = 0;                    % Zero the value of x_vl2 where the corresponding element is infinity

    x_isnan = isnan(x_vl2);              % Find elements of x_vl2 that are not of type numeric
    x_nan = find(x_isnan == 1);
    x_vl2(x_nan) = 0;                    % Set the value of elements in x_vl2 that are not of type numeric to zero
    
    % The corrected value is obtained according to the ipfp algorithm and saved
    x_val = ipfp(x_vl2, y_cur_vl, A, x_erro, y_erro, n_flow, n_link, n_loop, lb, []);
    IPFP(:, tm) = x_val;    
end 

close( h ) ;
