function [data_center,column_mean,row_mean,overall_mean,test_center,column_mean_test,row_mean_test]=doublecentering(data,test)

% Double-centering of distance/similarity/proximity matrices
% Input:
% data - training distance/similarity/proximity matrix
% test - test distance/similarity/proximity matrix
% Output:
% data_center - double-centered training distance/similarity/proximity matrix
% column_mean - training distance/similarity/proximity matrix column mean
% array
% array (with dimension suitable for double-centering the test distance/similarity/proximity matrix)
% row_mean - training distance/similarity/proximity matrix row mean array
% overall_mean - training distance/similarity/proximity matrix grand mean
% test_center - double-centered test distance/similarity/proximity matrix
% column_mean_test - training distance/similarity/proximity matrix column mean
% row_mean_test - test distance/similarity/proximity matrix row mean array

cm=ones(size(data,1),1)*mean(data);
rm=mean(data,2)*ones(1,size(data,2));
om=mean(data(:));

if nargin>1
    
    cmt=ones(size(test,1),1)*mean(data);
    rmt=mean(test,2)*ones(1,size(data,2));
    column_mean_test=cmt;
    row_mean_test=rmt;
    test_center=test-cmt-rmt+om;
    
end

data_center=data-cm-rm+om;
column_mean=cm;
row_mean=rm;
overall_mean=om;