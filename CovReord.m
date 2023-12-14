% C = [ "XX" "XY" "XX1" "XY1" "XX2" "XY2"; 
%     "XY" "YY" "YX1" "YY1" "YX2" "YY2";
%     "XX1" "YX1" "X1X1" "X1Y1" "X1X2" "X1Y2";
%     "XY1" "YY1" "X1Y1" "Y1Y1" "Y1X2" "Y1Y2";
%     "XX2" "YX2" "X1X2" "Y1X2" "X2X2" "X2Y2";
%     "XY2" "YY2" "X1Y2" "Y1Y2" "X2Y2" "Y2Y2";];

% B = CovReord(C, [2 2 2]);
% disp(B);


function C = CovReord(C, vs)
% reordering covariance matrix in the following order:
% X0 (first predictor), X1 (second predictor), S (target)

% indeces
odd = 1:2:sum(vs);
even = 2:2:sum(vs);
xx = 1:size(C,1)/2;
yy = size(C,1)/2+1:size(C,1);

Cov = C;
Cov2 = C;

% moving rows
Cov2(xx,:) = C(odd, :);
Cov2(yy,:) = C(even, :);

% moving columns
Cov(:, xx) = Cov2(:, odd);
Cov(:, yy) = Cov2(:, even);

% moving target at the end and reordering (joint future)
row_X = Cov(1,:);
row_Y = Cov(size(Cov,1)/2+1,:);

row_X = [row_X(2:length(row_X)/2) row_X(length(row_X)/2+2:end) ...
         row_X(1) row_X(length(row_X)/2+1)];

row_Y = [row_Y(2:length(row_Y)/2) row_Y(length(row_Y)/2+2:end) ...
         row_Y(1) row_Y(length(row_Y)/2+1)];

C = vertcat ( horzcat(Cov(2:size(Cov,1)/2,2:size(Cov,1)/2), ...
              Cov(2:size(Cov,1)/2,size(Cov,1)/2+2:end)), ...
              horzcat(Cov(size(Cov,1)/2+2:end,2:size(Cov,1)/2), ...
              Cov(size(Cov,1)/2+2:end,size(Cov,1)/2+2:end)));
            
C = horzcat(C, row_X(1:end-2)', row_Y(1:end-2)');
C = vertcat(C, row_X, row_Y);

end