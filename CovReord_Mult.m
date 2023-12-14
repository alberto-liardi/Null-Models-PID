% C = [ "XX" "XY" "XX1" "XY1" "XX2" "XY2"; 
%     "XY" "YY" "YX1" "YY1" "YX2" "YY2";
%     "XX1" "YX1" "X1X1" "X1Y1" "X1X2" "X1Y2";
%     "XY1" "YY1" "X1Y1" "Y1Y1" "Y1X2" "Y1Y2";
%     "XX2" "YX2" "X1X2" "Y1X2" "X2X2" "X2Y2";
%     "XY2" "YY2" "X1Y2" "Y1Y2" "X2Y2" "Y2Y2";];

% C = [ "XX" "XX'" "XY" "XY'" "XX1" "XX'1" "XY1" "XY'1"; 
%     "X'X" "X'X'" "X'Y" "X'Y'" "X'X1" "X'X'1" "X'Y1" "X'Y'1";
%     "YX" "YX'" "YY" "YY'" "YX1" "YX'1" "YY1" "YY'1";
%     "Y'X" "Y'X'" "Y'Y" "Y'Y'" "Y'X1" "Y'X'1" "Y'Y1" "Y'Y'1";
%     "X1X" "X1X'" "X1Y" "X1Y'" "X1X1" "X1X'1" "X1Y1" "X1Y'1"; 
%     "X'1X" "X'1X'" "X'1Y" "X'1Y'" "X'1X1" "X'1X'1" "X'1Y1" "X'1Y'1";
%     "Y1X" "Y1X'" "Y1Y" "Y1Y'" "Y1X1" "Y1X'1" "Y1Y1" "Y1Y'1";
%     "Y'1X" "Y'1X'" "Y'1Y" "Y'1Y'" "Y'1X1" "Y'1X'1" "Y'1Y1" "Y'1Y'1";];

% C = [ "XX" "XX'" "XX''" "XY" "XY'" "XX1" "XX'1" "XX''1" "XY1" "XY'1"; 
%     "X'X" "X'X'" "X'X''" "X'Y" "X'Y'" "X'X1" "X'X'1" "X'X''1" "X'Y1" "X'Y'1";
%     "X''X" "X''X'" "X''X''" "X''Y" "X''Y'" "X''X1" "X''X'1" "X''X''1" "X''Y1" "X''Y'1";
%     "YX" "YX'" "YX''" "YY" "YY'" "YX1" "YX'1" "YX''1" "YY1" "YY'1";
%     "Y'X" "Y'X'" "Y'X''" "Y'Y" "Y'Y'" "Y'X1" "Y'X'1" "Y'X''1" "Y'Y1" "Y'Y'1";
%     "X1X" "X1X'" "X1X''" "X1Y" "X1Y'" "X1X1" "X1X'1" "X1X''1" "X1Y1" "X1Y'1"; 
%     "X'1X" "X'1X'" "X'1X''" "X'1Y" "X'1Y'" "X'1X1" "X'1X'1" "X'1X''1" "X'1Y1" "X'1Y'1";
%     "X''1X" "X''1X'" "X''1X''" "X''1Y" "X''1Y'" "X''1X1" "X''1X'1" "X''1X''1" "X''1Y1" "X''1Y'1";
%     "Y1X" "Y1X'" "Y1X''" "Y1Y" "Y1Y'" "Y1X1" "Y1X'1" "Y1X''1" "Y1Y1" "Y1Y'1";
%     "Y'1X" "Y'1X'" "Y'1X''" "Y'1Y" "Y'1Y'" "Y'1X1" "Y'1X'1" "Y'1X''1" "Y'1Y1" "Y'1Y'1";];

% C
% C = magic(22);
% 
% p = 10;
% S1 = 1; 
% S2 = 1;
% 
% A = CovReord_Mult(C,S1,S2,p);
% disp(A);
% 

function C = CovReord_Mult(C,L1,L2,p)
% reordering covariance matrix in the following order:
% S1 (first source), S2 (second source), T (target)

% need to add some controls over dimensions...
assert( size(C,1) == size(C,2), "Covariance matrix is not a square matrix!");
assert( size(C,1) == (L1+L2)*(p+1), ...
        "Dimension of covariance matrix and sources-targets are not compatible");

vs = [L1*p L2*p L1+L2];
L = L1+L2;

% indeces
odd = 1:2:sum(vs);
even = 2:2:sum(vs);
xx = [1:L1*(p+1)];
yy = [L1*(p+1)+1:L1*(p+1)+L2*(p+1)];


% first I need to move S1 to the left, and S2 to the right
% target T is still at the beginning, will move it later

new_odd = zeros(L1,length(odd));
new_even = zeros(L2,length(even));

for j = 1:length(odd)
    new_odd(:,j) = (odd(j)-1)/2*L + [1:L1];
end
for j = 1:length(even)    
    new_even(:,j) = (even(j)-2)/2*L + L1 + [1:L2];
end

odd = reshape(new_odd, 1,size(new_odd,1)*size(new_odd,2));
even = reshape(new_even, 1,size(new_even,1)*size(new_even,2));

odd = odd(1:L1*(p+1));
even = even(1:L2*(p+1));

Cov = C;
Cov2 = C;
Cov3 = C; 

% moving rows
Cov2(xx,:) = C(odd, :);
Cov2(yy,:) = C(even, :);

% moving columns
Cov(:, xx) = Cov2(:, odd);
Cov(:, yy) = Cov2(:, even);

% disp("T1 - S1 - T2 - S2"); Cov

% now need to reorder sources among themselves 
% (have a format like S1_X S1_Y ... S1_Z, i.e. must be separated)

% Source 1
Cov3 = Cov;
l = 0:L1;
for j = 1:L1   
    
    pos = l(j)*(p+1)+1:1:(l(j)+1)*(p+1);
    index = [j:L1:L1*(p+1)];

    % moving rows
    Cov3(pos,:) = Cov(index,:);
end

% disp("move S1 rows"); Cov3

for j = 1:L1   
    
    pos = l(j)*(p+1)+1:1:(l(j)+1)*(p+1);
    index = [j:L1:L1*(p+1)];

    % moving columns
    Cov(:,pos) = Cov3(:, index);
end

Cov(:,L1*(p+1)+1:end) = Cov3(:, L1*(p+1)+1:end);

% disp("move S1 col"); Cov


% Source 2
Cov2 = Cov;
l = 0:L2;
for j = 1:L2   
    
    pos = L1*(p+1)+l(j)*(p+1)+1:1:(l(j)+1)*(p+1)+L1*(p+1);
    index = [L1*(p+1)+j:L2:L*(p+1)];
%problem here somehwere
    % moving rows
    Cov2(pos,:) = Cov(index,:);
end

% disp("move S2 rows"); Cov2

for j = 1:L2   
    
    pos = L1*(p+1)+l(j)*(p+1)+1:1:(l(j)+1)*(p+1)+L1*(p+1);
    index = [L1*(p+1)+j:L2:L*(p+1)];

    % moving columns
    Cov(:,pos) = Cov2(:, index);
end

Cov(:,1:L1*(p+1)) = Cov2(:,1:L1*(p+1));

% disp("move S2 cols"); Cov

C = Cov;

% moving target at the end and reordering (joint future)

T1 = 1:p+1:L1*(p+1);
T2 = L1*(p+1)+1:p+1:L*(p+1);

row_T1 = Cov(T1,:); % take T1 rows
row_T2 = Cov(T2,:); % take T2 rows

dim = L*(p+1);

% all targets
T = horzcat(T1, T2);

% avoiding targets
nT = setdiff([1:dim], T);

row_T1 = horzcat( row_T1(:,nT), row_T1(:,T));
row_T2 = horzcat( row_T2(:,nT), row_T2(:,T));

C = C(nT, nT);

C = horzcat(C, row_T1(:, 1:end-L)', row_T2(:, 1:end-L)');
C = vertcat(C, row_T1, row_T2);

end