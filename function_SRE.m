%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OMP algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% input: 
%       TM : measurements
%       OD : 
%
% output: 
%       SRE: Spatial relative errors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SRE = function_SRE(TM, real_od)
% Similar to the temporal relative error, find the spatial relative error of the prediction results
[r, c] = size(TM);
avgval_od = mean(real_od');
% Sorted by average of traffic size from smallest to largest
[smalltobigval, index] = sort(avgval_od);

res_realod = real_od(index, :);
res_x_TM = TM(index, :);
% Get the difference between the predicted result and the true traffic size matrix
res_TM = res_x_TM - res_realod;
% Relative error in space based on od stream numbering
for i = 1:r
    colnorm_realod(i) = norm(res_realod(i, :));
    colnorm_x_TM(i) = norm(res_TM(i, :));
end

errsp_x_TM = colnorm_x_TM ./ colnorm_realod;
% Eliminating irrational data from the data
% An exception is thrown when there are zeros in the divisor data
zero_x_TM = find(colnorm_x_TM == 0.0);
errsp_x_TM(zero_x_TM) = 0;

% return value
SRE = errsp_x_TM;
