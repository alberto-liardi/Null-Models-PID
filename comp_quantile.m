%% Computing the quantile of a given value in a distibution

% data = MxN matrix: each row is a distribution
% value = MxL matrix: the m-th row is compared with the m-th row of 'data'

function quantile = comp_quantile(data,value)

%     disp(size(data)); disp(size(value)); 
    assert(size(data,1)==size(value,1), "Dimensions not consistent!");
    
    quantile = zeros(size(value));

    for l = 1:size(value,2)
        nless = sum(data < value(:,l) - 1e-8,2);
        nequal = sum(abs(data-value(:,l))<1e-8,2);
        quantile(:,l) = (nless + 0.5.*nequal) / length(data);
    end

%     disp(size(quantile));
    
end