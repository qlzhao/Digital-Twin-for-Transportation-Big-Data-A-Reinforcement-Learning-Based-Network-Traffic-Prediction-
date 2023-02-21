function varargout = ipfp(in_x_vector, in_y_vector, in_a_matrix, in_err_x, in_err_y, n_x_val, n_y_val, n_loop, in_lb_vector, in_ub_vector)

% in_x_vector : TM(:, tm), A vector containing the size of all od flows at a given point in time
% in_y_vector : val_link(:, tm), Link load matrix, the load of all links at a given moment
% in_a_matrix : Routing Matrix, The rows represent all links of the network, and the columns represent all od flows, which, from a per-column perspective, signify which links each od flow has passed through
% in_err_x : x_erro(:, 1) = 0.001; OD flow traffic error bound
% in_err_y : y_erro(:, 1) = 0.001; link byte counts error bound
% n_x_val : n_flow = c_A; flow number
% n_y_val : n_link = r_A; link number
% n_loop : n_loop = 70; loop times
% in_lb_vector : zeros(c_A, 1); flow number
% in_ub_vector : []

% in_lb_vector indicates the minimum value of the traffic size, in the traffic size TM matrix, if the traffic size is less than the minimum value, then the value is assigned to the minimum value
if (length(in_lb_vector) > 0)
    x_islb = find(in_x_vector < in_lb_vector);
    in_x_vector(x_islb) = in_lb_vector(x_islb);
end

% Similarly, in_ub_vector represents the highest value of the traffic size, in the traffic size TM matrix, if the traffic size is greater than the maximum value, then assign this maximum value? The meaning of assigning 0 is unknown
if (length(in_ub_vector) > 0)
    x_isub = find(in_x_vector > in_ub_vector);
    in_x_vector(x_isub) = 0;
end

% Initialize boundary error variables
x_err = ones(n_x_val, 1); % flow number
y_err = ones(n_y_val, 1); % link number

% The num variable is the current maximum number of loop
num = 0;

% flag is the achievement flag of the condition (r_y1 == 0 & r_x1 == 0)
flag = 0;

% Enumerate the current traffic size vector
x_val = in_x_vector;

% Correction of the loop of the traffic size vector
while (1)
    % The initial traffic size vector is recorded, and the corrected traffic size vector is compared to this vector
    x_ini = x_val;
    
    % The second level for loop is link number of times
    for k = 1:n_y_val
        % in_a_matrix(k, :) is the passing of all od streams on the kth link, 1 if an od stream passes through the current link, 0 otherwise
        % l_vl is the sum of the traffic sizes of all od flows passing through the kth link
        l_vl = in_a_matrix(k, :) * x_val;
        % The value of the divisor of the current fraction cannot be 0
        if (l_vl == 0)
            continue;
        end
        
        % The real traffic size value of this link divided by the sum of the estimated traffic sizes of all od flows passing through this link
        % Get the current coefficient
        coef = in_y_vector(k) / l_vl;
        % Copy the current coefficients of coef into a vector of od flow numbers
        exco = repmat(coef, 1, n_x_val);
        % The traffic size is corrected by multiplying the traffic size value by the current correction factor if the od flow passes through the current link, otherwise the current traffic size value is not corrected
        exrs = in_a_matrix(k, :) .* exco + (1 - in_a_matrix(k, :));
        % The traffic size matrix is multiplied by the current coefficient correction matrix to obtain the new traffic size matrix
        x_val = x_val .* exrs';
    end
    
    % The following is a sign of early exit of the loop
    % Calculate the current link traffic size error, the core of this algorithm is to make the sum of the traffic size of all od flows passing through a link approximates the current link traffic size value
    y_err = in_a_matrix * x_val - in_y_vector;
    % Calculate whether the magnitude of correction of the traffic size vector is too small
    x_err = x_val - x_ini;
    y = abs(y_err);
    x = abs(x_err);
    
    % Determine if there is data greater than the error limit, if not then flag is set to 1 and exit the loop early
    ry = y > in_err_y;
    rx = x > in_err_x;
    y1 = find(ry == 1);
    x1 = find(rx == 1);

    [r_y1, c_y1] = size(y1);
    [r_x1, c_x1] = size(x1);
    if (r_y1 == 0 & r_x1 == 0)
        flag = 1;
        break;
    end
    
    % The number of cycles, equal to the maximum number of loops and then exit
    num = num + 1;
    if (num == n_loop)
        break;
    end
end

% nargout indicates the number of output arguments, and the return value is obtained based on the number of output arguments
switch nargout
    case 1
        varargout{1} = x_val;
    case 2
        varargout{1} = x_val;
        varargout{2} = num;
    case 3
        varargout{1} = x_val;
        % Number of iterations
        varargout{2} = num;
        % Markers for whether an iteration is completed early
        varargout{3} = flag;
    case 4
        varargout{1} = x_val;
        varargout{2} = num;
        varargout{3} = flag;
        varargout{4} = x_err;
    case 5
        varargout{1} = x_val;
        varargout{2} = num;
        varargout{3} = flag;
        % OD flow traffic error bound
        varargout{4} = x_err;
        % link byte counts error bound
        varargout{5} = y_err;
end
