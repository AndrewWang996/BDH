fC2R = @(x) [real(x) imag(x)];
fR2C = @(x) complex(x(:,1), x(:,2));

%%
% shapename = '4.Trl5_Camera01';
% shapename = 'Monkey_square2';
% shapename = 'red_dragon_square_t2';
% imgfilepath = [cd '\data\' shapename '.png'];
% cagefilepath = [cd '\data\' shapename '_cage.obj'];

if exist('datadir', 'var')~=1 
    datadir = [cd '\data\red dragon'];
    warning('data path not set, default to %s', datadir);
end

if exist('numMeshVertex', 'var')~=1 
    numMeshVertex = 10000;
    warning('numMeshVertex not set, default to %d!', numMeshVertex);
end

imgfilepath  = [datadir '\image.png'];
cagefilepath = [datadir '\cage.obj'];
datafile = [datadir 'data.mat'];

%%
% cage = fR2C( readObj( cagefilepath ) );
[cx, cf] = readObj( cagefilepath );
cage = fR2C(cx(cf,:));
if abs(cage(1)-cage(end))<1e-3, cage = cage(1:end-1); end
if signedpolyarea(cage)<0; cage = cage(end:-1:1); end

P2Psrc = zeros(0,1); P2Pdst = zeros(0,1);

if exist(datafile, 'file') == 2
    load(datafile);
    hasDataLoaded = true;
end


% P2Psrc=P2Psrc([1:5 8:end]);
% P2Pdst(4:7)=P2Psrc(4:7);

% cage_offset = 0.1;

% d = distance2polygon(real(P2Psrc), imag(P2Psrc), real(cage), imag(cage));

if exist([datadir 's.obj'], 'file') == 2
    fprintf('loading mesh from s.obj\n');
    [X, T] = readObj([datadir 's.obj']);
else
    [X, T] = cdt(fC2R(cage), [], numMeshVertex, false); % 50000 as in maya, replace the following line, 'cause the P2Psrc may fall on the cage, case problematic triangulation (very small triangles)
end

% [X, T] = cdt(fC2R(cage), fC2R(P2Psrc), numMeshVertex, false); % 50000 as in maya
% [X, T] = genRegularMesh(fC2R(cage), numMeshVertex);
% [X, T] = cdt(fC2R(cage), fC2R(P2Psrc), 50000, 0); % 50000 as in maya
X = fR2C(X);

%% load handles
% P2PVtxIds = numel(cage) + (1:numel(P2Psrc));
% P2PVtxIds = triangulation(T, fC2R(X)).nearestNeighbor( fC2R(P2Psrc) );
P2PVtxIds = zeros(numel(P2Psrc), 1);
for i=1:numel(P2Psrc),   [~, P2PVtxIds(i)] = min( abs(X-P2Psrc(i)) );  end


% [X, T] = cdt(fC2R(cage), [], numMeshVertex, false); % 50000 as in maya
% X = fR2C(X);
% 
% %% load handles 
% P2PVtxIds = zeros(numel(P2Psrc), 1);
% for i=1:numel(P2Psrc)
%     [~, P2PVtxIds(i)] = min( abs(X-P2Psrc(i)) );
% end
P2PCurrentPositions = P2Pdst;


%% texture
% img = imread(imgfilepath, 'BackgroundColor', [1 1 1 1]);
img = imread(imgfilepath);
[w, h, ~] = size(img);
uv = fR2C([real(X)/h imag(X)/w])*100 + complex(0.5, 0.5);


fDrawPoly = @(x, c) plot(real(x([1:end 1])), imag(x([1:end 1])), c);


%%
% fPlot = @(x, c) plot(real(x), imag(x), c);
% figuredocked;  fPlot(v, 'r-');
% hold on;  fPlot(offsetcage, 'b-');
% hm = drawmesh(t,x);

% !start glvu_bdh.ex_

%% for Key Frame interpolation
if exist([datadir 'PhiPsyKF.mat'], 'file') == 2
    load([datadir 'PhiPsyKF']);
    fprintf('%d frames are loaded\n', size(PhiPsyKF,2)/2);
    
    if exist('anchorsForInterp', 'var') == 1 && ~isempty(anchorsForInterp)
%         interpAnchID = triangulation(T, fC2R(X)).nearestNeighbor( fC2R(anchorsForInterp) );
        interpAnchID = zeros(numel(anchorsForInterp), 1);
        for i=1:numel(anchorsForInterp),   [~, interpAnchID(i)] = min( abs(X-anchorsForInterp(i)) );  end
    end
end

ikeyframe = 2;
BDHIAlpha = 1;

offsetCage = v;

if ~hasGPUComputing, warning('NO CUDA capable GPU is present, the computation will fall back to CPU!'); end
