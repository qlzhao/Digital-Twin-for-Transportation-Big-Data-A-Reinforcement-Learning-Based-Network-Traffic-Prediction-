%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OMP algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% input: 
%       TM : measurements
%       OD : 
%
%output: 
%       TRE: Temporal relative errors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TRE = function_TRE(TM, real_od)
% Rows represent od flows, columns represent time scales
[r, c] = size(real_od);
% Get the difference between the predicted flow and the true flow size
res_diffod_TM = TM - real_od;

% Solve for each moment, with the numerator being the two-parameter number of the difference and the denominator being the two-parameter number of the true value
for i = 1:c
    % The norm function is used to find the norm of the matrix and defaults to the second norm.
    rownorm_realod(i) = norm(real_od( :, i));
    rownorm_diffod_TM(i) = norm(res_diffod_TM(:, i));
end

errtm_vector_TM = rownorm_diffod_TM ./ rownorm_realod;
TRE = errtm_vector_TM;
