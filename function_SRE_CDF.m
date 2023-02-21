% Spatial Relative Error CDF
function [errsp_vector_TM_1, errsp_vector_TM_cdf_1] = function_SRE_CDF(errsp_TM)
% Sorted by spatial relative error from smallest to largest
[errsp_vector_TM, ind_TM] = sort(errsp_TM);
% Intermediate variables
i = 1; k = 1; sum = 0;
ll = length(errsp_vector_TM);

while( 1 )
    % while loop method, traversing from the first index to the last index
    if (i > ll )
        break;
    end
    
    % Iteration of error values from smallest to largest
    val = errsp_vector_TM(i);
    a(k) = val;
    % Count the number of error with a val value
    num1 = find(errsp_vector_TM == val);
    num2 = length(num1);
    % The scale of the coordinates of the y-axis is constantly accumulating
    sum = sum + num2 / ll;
    b(k) = sum;
    i = i + num2;
    k = k + 1;
end

% return value, The x and y axes of the probability density function
errsp_vector_TM_1 = a;
errsp_vector_TM_cdf_1 = b;
