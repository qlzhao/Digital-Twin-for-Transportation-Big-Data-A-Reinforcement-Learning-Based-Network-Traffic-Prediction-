% Temporal Relative Error CDF
function [errtm_vector_val_TM_1, errtm_vector_TM_cdf_1] = function_TRE_CDF(errtm_vector_TM)
% Sorted by Temporal relative error from smallest to largest
[errtm_vector_val_TM, ind_haar] = sort(errtm_vector_TM);
% Obtain the number of rows and columns of the Temporal relative error matrix
[r_errtm_vector_TM, c_errtm_vector_TM] = size(errtm_vector_TM);
% a records the x-axis of the probability density function and b records the y-axis of the probability density function
a = [];
b = [];
% Intermediate variables
i = 1; k = 1; sum = 0;
ll = length(errtm_vector_val_TM);

while( 1 )
    % while loop method, traversing from the first index to the last index
    if (i > ll)
        break;
    end
    % Iteration of error values from smallest to largest
    val = errtm_vector_val_TM(i);
    a(k) = val;
    % Count the number of error with a val value
    num1 = find(errtm_vector_val_TM == val);
    num2 = length(num1);
    sum = sum + num2 / ll;
    % The scale of the coordinates of the y-axis is constantly accumulating
    b(k) = sum;
    i = i + num2;
    k = k + 1;
end

% return value
errtm_vector_val_TM_1 = a;
errtm_vector_TM_cdf_1 = b;
