%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This script is a rough version for sampling PF points from DC-DTLZ
% problems. Since the constraints in DC-DTLZ problems are not only about
% objective vectors, it is a bit difficult to get a closed form to derive
% the uniformly distributed PF samples.
%
% Author: Ke Li
% E-mail: k.li@exeter.ac.uk || keli.genius@gmail.com
% University of Exeter
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% parameter settings

close all
clear
clc

objDim = 3;

a = 5;
b = 0.5;   % b = 0.95 for DC1 and b = 0.5 for DC3

dc_id  = 3; % constraint type
p_id   = 2; % p_id = 1 for DTLZ1 and 2 for DTLZ2-4
dim    = objDim;
lowend = zeros(1, dim);
span   = ones(1, dim) - lowend;

sample_size    = 10000; % this number needs to be adapted according to your requirement and the number of objectives. HINT: obviously, in a high-D space, the number of samples should be large enough
random_samples = lowend(ones(1, sample_size), :) + span(ones(1, sample_size), :) .* lhsdesign(sample_size, dim, 'criterion', 'maximin', 'iterations', 1000);
y_samples      = zeros(size(random_samples, 1), objDim);

%% sample PF points from non-constraint space
if p_id == 1    % DTLZ1
    t   = 1;
    sum = t;
    for i = 1 : objDim - 1
        sum = sum .* random_samples(:, i);
    end
    y_samples(:, 1) = 0.5 .* sum;

    for i = 2 : objDim
        sum = t;
        for j = 1 : objDim - i
            sum = sum .* random_samples(:, j);
        end
        sum             = sum .* (1.0 - random_samples(:, objDim - i + 1));
        y_samples(:, i) = 0.5 .* sum;
    end
elseif p_id == 2    % DTLZ2-4
    t   = 1;
    sum = t;
    for i = 1 : objDim - 1
        sum = sum .* cos(random_samples(:, i) .* pi ./ 2.0);
    end
    y_samples(:, 1) = sum;

    for i = 2 : objDim
        sum = t;
        for j = 1 : objDim - i
            sum = sum .* cos(random_samples(:, j) .* pi ./ 2.0);
        end
        sum = sum .* sin(random_samples(:, objDim - i + 1) .* pi ./ 2.0);
        y_samples(:, i) = sum;
    end
end

%% remove infeasible solutions
if dc_id == 1
    con = cos(a * pi * random_samples(:, 1)) - b;
    idx = find(con > 0);
    PF  = y_samples(idx, :);
elseif dc_id == 2
    PF = y_samples;
elseif dc_id == 3
    diff = cos(a * pi * random_samples) - b;
    con  = min(diff, [], 2);
    idx  = find(con > 0);
    PF   = y_samples(idx, :);
else
    error('Undefined constraint type!')
end

%% plot results
if objDim == 2
    plot(PF(:, 1), PF(:, 2), 'o');
elseif objDim == 3
    plot3(PF(:, 1), PF(:, 2), PF(:, 3), 'o');
else
    parallelcoords(PF);
end
save('DC3_DTLZ3.pf','PF','-ascii');