fprintf('\n\n');
fprintf('we are running p2p_harmonic, bro');

active_set_is_on = true; %should we use the active set approach or not
active_set_visualization = false; %visualization of active set on/off
maxBisectionIterations = 10;

use_Cauchy_argument_principle = false; %using true may slows down the validation - not needed in general. use false

%preprocessing - done once
if (needsPreprocessing)
%if(~exist('E', 'var'))

    if exist('Phi', 'var')~=1 || numel(Phi) ~= numel(v)
        Phi = v;
        Psy = zeros(size(Phi));
    end
    Phi_start = Phi;
    

	%preprocess for Modulus of Continuity
	[L, indexOfMaxL] = computeLipschitzConstantOfDerivativeOfCauchy(v, denseEvaluationSamples);
    if hasGPUComputing
        L = gpuArray(L);
        indexOfMaxL = gpuArray(indexOfMaxL);

        v = gpuArray(v);
        
        denseEvaluationSamples = gpuArray(denseEvaluationSamples);
    end

    DerivativeOfCauchyCoordinatesAtFixedSamples = gather(derivativesOfCauchyCoord(v, fixedSamples.')); %copmute on gpu and transfer to cpu
    DerivativeOfCauchyCoordinatesAtEnergySamples = gather(derivativesOfCauchyCoord(v, energySamples.')); %copmute on gpu and transfer to cpu
    DerivativeOfCauchyCoordinatesAtActiveSetPoolSamples = gather(derivativesOfCauchyCoord(v, activeSetPoolSamples.')); %copmute on gpu and transfer to cpu
    DerivativeOfCauchyCoordinatesAtDenseEvaluationSamples = derivativesOfCauchyCoord(v, denseEvaluationSamples.'); %no gather - keep on the gpu
    
    if exist('P2Ppositions', 'var')==1
        CauchyCoordinatesAtP2Phandles = cauchyCoordinates(v, P2Ppositions.'); %compute on the cpu since the number of P2Ps is very small
    end

    fillDistanceSegments = abs(circshift(denseEvaluationSamples, -1) - denseEvaluationSamples)/2; %fill distance for each segment
    
    %preprocess for ARAP energy
    Q = DerivativeOfCauchyCoordinatesAtEnergySamples'*DerivativeOfCauchyCoordinatesAtEnergySamples;
    %sqrtQ = sqrtm(Q);
    
    mydiag = @(x) sparse(1:numel(x), 1:numel(x), x);
    [E, vec] = eig(Q);
    vec = diag(vec); vec = vec(2:end);
    ARAP_q = mydiag(vec.^0.5)*E(:,2:end)';
    clear Q;
    
    needsPreprocessing = false;

    forceConformalMode = false;
end

total_time = tic;

if(active_set_is_on)
    
    if(~exist('activeSetSigma1', 'var'))
        activeSetSigma1 = [];
    end
    if(~exist('activeSetSigma2', 'var'))
        activeSetSigma2 = [];
    end
    if(~exist('activeSet_k', 'var'))
        activeSet_k = [];
    end
    
    abs_fz = abs(DerivativeOfCauchyCoordinatesAtActiveSetPoolSamples*Phi);
    abs_fzbar = abs(DerivativeOfCauchyCoordinatesAtActiveSetPoolSamples*Psy);
    
    sigma1 = abs_fz + abs_fzbar;
    sigma2 = abs_fz - abs_fzbar;
    k = abs_fzbar ./ abs_fz;

    upperThresholdSigma1 = 0.95*sigma1_upper_bound;
    lowerThresholdSigma1 = 0.945*sigma1_upper_bound;

    lowerThresholdSigma2 = 1.15*sigma2_lower_bound;
    upperThresholdSigma2 = 1.2*sigma2_lower_bound;

    upperThreshold_k = 0.95*k_upper_bound;
    lowerThreshold_k = 0.945*k_upper_bound;

    warning('off', 'signal:findpeaks:largeMinPeakHeight');
    
    [~, indicesToBeAddedSigma1] = findpeaks(sigma1, 'MinPeakHeight', upperThresholdSigma1);
    indicesToBeAddedSigma1 = union(indicesToBeAddedSigma1, find(sigma1 > sigma1_upper_bound));
    [~, indicesToBeAddedSigma2] = findpeaks(-sigma2, 'MinPeakHeight', -lowerThresholdSigma2); %we use minus to find local minimum
    indicesToBeAddedSigma2 = union(indicesToBeAddedSigma2, find(sigma2 < sigma2_lower_bound));
    [~, indicesToBeAdded_k] = findpeaks(k, 'MinPeakHeight', upperThreshold_k);
    indicesToBeAdded_k = union(indicesToBeAdded_k, find(k > k_upper_bound));

        
    activeSetSigma1 = union(activeSetSigma1, indicesToBeAddedSigma1);
    indicesToBeRemovedSigma1 = activeSetSigma1(sigma1(activeSetSigma1) < lowerThresholdSigma1);
    activeSetSigma1 = setdiff(activeSetSigma1, indicesToBeRemovedSigma1);
    
    activeSetSigma2 = union(activeSetSigma2, indicesToBeAddedSigma2);
    indicesToBeRemovedSigma2 = activeSetSigma2(sigma2(activeSetSigma2) > upperThresholdSigma2);
    activeSetSigma2 = setdiff(activeSetSigma2, indicesToBeRemovedSigma2);
    
    activeSet_k = union(activeSet_k, indicesToBeAdded_k);
    indicesToBeRemoved_k = activeSet_k(k(activeSet_k) < lowerThreshold_k);
    activeSet_k = setdiff(activeSet_k, indicesToBeRemoved_k);

    
    numActiveSetPool = size(DerivativeOfCauchyCoordinatesAtActiveSetPoolSamples, 1);

    if(active_set_visualization)
        %figure;

        subplot(2, 1, 1);
        plot(sigma1, 'Color', 'black');
        axis([0 numActiveSetPool sigma2_lower_bound-0.5 sigma1_upper_bound+0.5]);
        title(['{' '\color{black} \fontsize{10}' '# sigma1 samples:' '\color{cyan} \fontsize{14}'  num2str(size(activeSetSigma1, 1)) '\color{black} \fontsize{10}' ', # sigma2 samples:' '\color{magenta} \fontsize{14}'  num2str(size(activeSetSigma2, 1)) '}'])
        xlabel('samples');
		ylabel('sigma1 & sigma2');
        hold on;
        plot(activeSetSigma1, sigma1(activeSetSigma1), '.', 'MarkerEdgeColor', 'cyan');
        plot([0 numActiveSetPool], [sigma1_upper_bound sigma1_upper_bound], 'Color', 'red', 'LineStyle', ':');

        plot(sigma2, 'Color', 'blue');
        plot(activeSetSigma2, sigma2(activeSetSigma2), '.', 'MarkerEdgeColor', 'magenta');
        plot([0 numActiveSetPool], [sigma2_lower_bound sigma2_lower_bound], 'Color', 'red', 'LineStyle', ':');
        hold off;

        subplot(2, 1, 2);
        plot(k);
		axis([0 numActiveSetPool -0.1 k_upper_bound+0.1]);
        title(['{# k samples:' '\color{green} \fontsize{14}'  num2str(size(activeSet_k, 1)) '}'])
        xlabel('samples');
		ylabel('k');
        hold on;
        plot(activeSet_k, k(activeSet_k), '.', 'MarkerEdgeColor', 'green');
        plot([0 numActiveSetPool], [k_upper_bound k_upper_bound], 'Color', 'red', 'LineStyle', ':');
        plot([0 numActiveSetPool], [0 0], 'Color', 'red', 'LineStyle', ':');
        hold off;
    end
    
else
    activeSetSigma1 = [];
    activeSetSigma2 = [];
    activeSet_k = [];
end

DerivativeOfCauchyCoordinatesAtActiveSamplesSigma1 = DerivativeOfCauchyCoordinatesAtActiveSetPoolSamples(activeSetSigma1, 1:end);
DerivativeOfCauchyCoordinatesAtActiveSamplesSigma2 = DerivativeOfCauchyCoordinatesAtActiveSetPoolSamples(activeSetSigma2, 1:end);
DerivativeOfCauchyCoordinatesAtActiveSamples_k = DerivativeOfCauchyCoordinatesAtActiveSetPoolSamples(activeSet_k, 1:end);

frames_sigma2 = calc_frames(DerivativeOfCauchyCoordinatesAtActiveSamplesSigma2, Phi);
frames_k = calc_frames(DerivativeOfCauchyCoordinatesAtActiveSamples_k, Phi);


numVirtualVertices = size(CauchyCoordinatesAtP2Phandles, 2);
numFixedSamples = size(DerivativeOfCauchyCoordinatesAtFixedSamples, 1);
numEnergySamples = size(DerivativeOfCauchyCoordinatesAtEnergySamples, 1);
numDenseEvaluationSamples = size(DerivativeOfCauchyCoordinatesAtDenseEvaluationSamples, 1);



frames_fixed = calc_frames(DerivativeOfCauchyCoordinatesAtFixedSamples, Phi);
%frames_fixed = ones(size(DerivativeOfCauchyCoordinatesAtFixedSamples, 1), 1); %reset the frames

frames_energy = calc_frames(DerivativeOfCauchyCoordinatesAtEnergySamples, Phi);

arap_frames_vector = transpose(frames_energy)*DerivativeOfCauchyCoordinatesAtEnergySamples; %note the use of transpose rather than '
% g = arap_frames_vector*E(:,2:end)*mydiag(vec.^-0.5);
ARAP_g = (arap_frames_vector*E(:,2:end)).*reshape(vec.^-0.5, 1, []);


optimizationTimeOnly = tic;

switch solver_type
    case 'CVX' 
        [solverStatus, Energy_total, E_ARAP, E_POSITIONAL, phi, psy] = cvx_p2p_harmonic( ...
            CauchyCoordinatesAtP2Phandles, DerivativeOfCauchyCoordinatesAtFixedSamples, ...
            DerivativeOfCauchyCoordinatesAtActiveSamplesSigma1, DerivativeOfCauchyCoordinatesAtActiveSamplesSigma2, DerivativeOfCauchyCoordinatesAtActiveSamples_k, ...
            P2PCurrentPositions, ARAP_g, ARAP_q, frames_fixed, frames_sigma2, frames_k, ...
            p2p_weight, sigma2_lower_bound, sigma1_upper_bound, k_upper_bound, ...
            numVirtualVertices, numFixedSamples, numEnergySamples, ...
            no_output, forceConformalMode);

    case 'Direct Mosek'
        if forceConformalMode
           error('Conformal mode is currently supported only in CVX. Switch to CVX or try to reduce bound on k to a small number in order to approximate conformal');
        end

        [solverStatus, Energy_total, E_ARAP, E_POSITIONAL, phi, psy] = mosek_p2p_harmonic( ...
            CauchyCoordinatesAtP2Phandles, DerivativeOfCauchyCoordinatesAtFixedSamples, ...
            DerivativeOfCauchyCoordinatesAtActiveSamplesSigma1, DerivativeOfCauchyCoordinatesAtActiveSamplesSigma2, DerivativeOfCauchyCoordinatesAtActiveSamples_k, ...
            P2PCurrentPositions, ARAP_g, ARAP_q, frames_fixed, frames_sigma2, frames_k, ...
            p2p_weight, sigma2_lower_bound, sigma1_upper_bound, k_upper_bound, ...
            numVirtualVertices, numFixedSamples, numEnergySamples, ...
            no_output);
        
    case 'Least Square Conformal'
        
        lambda = 1e-2;
        A = [CauchyCoordinatesAtP2Phandles; lambda*SODerivativeOfCauchyCoordinatesAtFixedSamples]; b = [P2PCurrentPositions; zeros(numFixedSamples,1)];
        phi = A\b;
        psy = phi*0;
        E_ARAP = norm( A*phi-b );
        E_POSITIONAL = 0;
        Energy_total = E_ARAP + E_POSITIONAL;
        solverStatus = 'solved';

    otherwise
        warning('Unexpected solver type.');
end

fprintf('Optimization time only time:%.4f\n', toc(optimizationTimeOnly));

validationTime = tic;

if exist('binarySearchValidMap', 'var')==1 && binarySearchValidMap
    [Phi, Psy, t] = validateMapBounds(L, indexOfMaxL, v, fillDistanceSegments, Phi, phi, Psy, psy, DerivativeOfCauchyCoordinatesAtDenseEvaluationSamples, sigma2_lower_bound, sigma1_upper_bound, k_upper_bound, maxBisectionIterations, unrefinedInwardCage, use_Cauchy_argument_principle);

    fprintf('Validation time:%.4f\n', toc(validationTime));
else
    t = 1; Phi = phi; Psy = psy;
end

DeformedP2PhandlePositions = CauchyCoordinatesAtP2Phandles*Phi + conj(CauchyCoordinatesAtP2Phandles*Psy);


fprintf('Solver status: %s\n', solverStatus);
fprintf('E_ARAP: %.8f, E_POS: %.8f, E: %.8f\n', E_ARAP, E_POSITIONAL, Energy_total);
fprintf('Total script time:%.4f\n', toc(total_time));
fprintf('Fixed_samples:%d, Sigma1_samples:%d, Sigma2_samples:%d, k_samples:%d\n', numFixedSamples, size(activeSetSigma1, 1), size(activeSetSigma2, 1), size(activeSet_k, 1));
fprintf('Vertices:%d, Energy_samples:%d, Evaluation_samples:%d\n', numVirtualVertices, numEnergySamples, numDenseEvaluationSamples);

mink = min( abs( (D*psy)./(D*phi) ) );
fprintf('min k: %f\n', mink);
