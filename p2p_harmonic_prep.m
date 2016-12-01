fprintf('\n Running p2p_harmonic_prep');
%% set parameters

if exist('cage_offset', 'var')~=1
    cage_offset = 2e-2;
end

% cage_offset = 5e-2;

if exist('numVirtualVertices', 'var')~=1
    fprintf('setting default parameters for p2p-harmonic deformation\n');
    numVirtualVertices = 1;
    numFixedSamples = 1;
    numEnergySamples = 10000;
    numDenseEvaluationSamples = 15000;
    numActiveSetPoolSamples = 2000;
    solver_type = 'Direct Mosek';
    % solver_type = 'CVX';

    no_output = true;
    p2p_weight = 100;
    sigma2_lower_bound = 0.35;
    sigma1_upper_bound = 5; 
    k_upper_bound = 0.8;
    binarySearchValidMap = true;
end



%% compute cauchy coordinates and its derivatives for samples
offsetCage = subdivPolyMat(cage, numVirtualVertices-1)*polygonOffset(cage, -cage_offset, false);

fSampleOnCage = @(n) subdivPolyMat(cage, n)*cage;

energySamples = fSampleOnCage(numEnergySamples);
fixedSamples = fSampleOnCage(numFixedSamples-1);
denseEvaluationSamples = fSampleOnCage(numDenseEvaluationSamples);
activeSetPoolSamples = fSampleOnCage(numActiveSetPoolSamples);

% DerivativeOfCauchyCoordinatesAtFixedSamples = regularCauchyCoordDerivative(offsetCage, fixedSamples);
% DerivativeOfCauchyCoordinatesAtEnergySamples = regularCauchyCoordDerivative(offsetCage, energySamples);
% DerivativeOfCauchyCoordinatesAtDenseEvaluationSamples = regularCauchyCoordDerivative(offsetCage, denseEvaluationSamples);
% DerivativeOfCauchyCoordinatesAtActiveSetPoolSamples = regularCauchyCoordDerivative(offsetCage, activeSetPoolSamples);

% [~, ~, SODerivativeOfCauchyCoordinatesAtFixedSamples] = regularCauchyCoord(offsetCage, fixedSamples);
[~, SODerivativeOfCauchyCoordinatesAtFixedSamples] = derivativesOfCauchyCoord(offsetCage, fixedSamples.');



%%
v = offsetCage;
unrefinedInwardCage = cage;

% compute D for interpolation
% [C, D] = regularCauchyCoord(offsetCage, X);
C = cauchyCoordinates(offsetCage, X.');
D = derivativesOfCauchyCoord(offsetCage, X.');

internalPoints = X;

clear E;

numFixedSamples = numel(fixedSamples);
numVirtualVertices = numel(offsetCage);

needsPreprocessing = true;

%%
if ~exist('Phi', 'var') || numel(Phi)~=numVirtualVertices
%     Phi = offsetCage;
%     Psy = Phi*0;
end


%%
myinpoly = @(x, y) inpolygon(real(x), imag(x), real(y), imag(y));
assert( signedpolyarea( fC2R(offsetCage) ) > 0, 'cage vertex order reversed for proper Cauchy coordiates computation!');
if ~isempty(selfintersect(real(offsetCage), imag(offsetCage))), warning('source polygon has self-interesction'); end
% assert( all(myinpoly(w, offsetCage)), 'Boundary samples should be on offset inside the cage' );
% assert( isempty(selfintersect(real(w), imag(w))), 'offset polygon has self-interesction');

% fDrawPoly = @(x) drawmesh([1:size(x,1); 2:size(x,1) 1; 1:size(x,1)]', x); % can be colorcoded
