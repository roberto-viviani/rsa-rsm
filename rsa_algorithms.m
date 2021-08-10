function rsa = rsa_algorithms
% The namespace of the algorithms used to compute the rsa.
%
% 2020-21 Roberto Viviani - December 2020
% Institute of Psychology, University of Innsbruck
% Last modified: August 1st 2021

    %A factory function to create the callback that computes brain pattern
    %similarities and their correspondence to stimuli similarities
    rsa.corrfun = @create_corrfun_;
    %The algorithms to compose similarity matrices from data features
    rsa.create_rsm = @create_rsm_;
end

function corrf = create_corrfun_(covv, pcovv, trmsk, args)
    %the callbacks computing the correlation (closures). These callbacks
    %are a function of the brain searchlight signal. covv are the
    %offdiagonal elements of the stimuli similarity, pcovv the offdiagonal
    %elements of the confound similarities, trmsk selects the offdiagonal
    %elements from a matrix (taking into account an off-diagonal offset if 
    %required).

    ccovv = size(covv, 2);
    method = args.Method;
    defs = rsa_defaults();
    if args.maxslsize < 27,  minsize = args.maxslsize;
    else,                    minsize = defs.minsearchlightsize;
    end
    utils = rsa_utils();
    cov2cor = utils.cov2cor;
    
    %the callbacks.
    function vec = f(a, b)   %simple correlation with sscp
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        cvox = size(a, 2);
        if cvox < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        BB = a * a' ./ cvox;
        vec = corr(BB(trmsk), covv, 'type', method);
    end
    function vec = f_cov(a, b) %simple correlation with mean-adjusted sscp
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        cvox = size(a, 2);
        if cvox < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        a = bsxfun(@minus, a, mean(a, 2));
        BB = a * a' ./ cvox;
        vec = corr(BB(trmsk), covv, 'type', method);
    end
    function vec = f_cor(a, b) %simple correlation with correlation
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        if size(a, 2) < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        %the same as BB = corr(a');
        a = bsxfun(@minus, a, mean(a, 2));
        BB = cov2cor(a * a');
        vec = corr(BB(trmsk), covv, 'type', method);
    end

    function vec = pf(a, b)     %partial correlation with sscp
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        cvox = size(a, 2);
        if cvox < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        BB = a * a' ./ cvox;
        vec = partialcorr(BB(trmsk), covv, pcovv, 'type', method);
    end
    function vec = pf_cov(a, b) %partial correl. with mean-adjusted sscp
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        cvox = size(a, 2);
        if cvox < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        a = bsxfun(@minus, a, mean(a, 2));
        BB = a * a' ./ cvox;
        vec = partialcorr(BB(trmsk), covv, pcovv, 'type', method);
    end
    function vec = pf_cor(a, b) %partial correlation with correlation
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        if size(a, 2) < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        %the same as BB = corr(a');
        a = bsxfun(@minus, a, mean(a, 2));
        BB = cov2cor(a * a');
        vec = partialcorr(BB(trmsk), covv, pcovv, 'type', method);
    end

    prm = []; rcvect = 0;        %regression
    if strcmpi(method, 'regression')
        nrows = size(covv, 1);
        prm = pinv((eye(nrows) - ones(nrows)./nrows) * [covv, pcovv]);
        rcvect = size(prm, 1);
    end
    function vec = lf(a, b)
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        if size(a, 2) < minsize
            vec = NaN(rcvect, 1);
            return;
        end
        a = bsxfun(@minus, a, mean(a, 2));
        BB = cov2cor(a * a');
        BB = BB(trmsk);
        vec = prm * BB;
    end

    function vec = sf(a, b)     %number of voxels in searchlight
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        cvox = size(a, 2);
        if cvox < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        vec = ones(ccovv,1) .* cvox;
    end

    %other distance functions, working through dissimilarities.
    distfname = args.BrainMapType;
    dcovv = 1 - bsxfun(@rdivide, covv, max(covv)); 
    function vec = f_distf(a, b)     
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        if size(a, 2) < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        vec = pdist(a, distfname)';
        vec = corr(vec, dcovv, 'type', method);
    end

    dpcovv = 1 - bsxfun(@rdivide, pcovv, max(pcovv)); 
    function vec = f_distf_cov(a, b)    
        idx = any(isnan(a));
        if any(idx), a(:,idx) = []; end
        if size(a, 2) < minsize
            vec = NaN(ccovv, 1);
            return;
        end
        vec = pdist(a, distfname)';
        vec = partialcorr(vec, dcovv, dpcovv, 'type', method);
    end

    %callback selection
    distfs = {'euclidean', 'squaredeuclidean', 'seuclidean', 'cosine', ...
        'correlation', 'spearman', 'jaccard', 'mahalanobis'};
    if strcmpi(method, 'regression')
        corrf = @lf;
    elseif strcmpi(method, 'slsize')
        corrf = @sf;
    else
        if isempty(pcovv)
            switch args.BrainMapType
                case 'sscp'
                    corrf = @f;
                case 'cov'
                    corrf = @f_cov;
                case 'cor'
                    corrf = @f_cor;
                case distfs
                    corrf = @f_distf;
                otherwise
                    error(['Invalid similarity spec: ', args.BrainMapType]);
            end
        else
            switch args.BrainMapType
                case 'sscp'
                    corrf = @pf;
                case 'cov'
                    corrf = @pf_cov;
                case 'cor'
                    corrf = @pf_cor;
                case distfs
                    corrf = @f_distf_cov;
                otherwise
                    error(['Invalid similarity spec: ', args.BrainMapType]);
            end
        end
    end

end

function rsm = create_rsm_(data, distf)
% Create a representational similarity matrix. data can contain
% categorical labels (cell arrays of strings or int8, int16 types etc.) or
% quantitative features (float data). In this latter case, a distance
% function may be specified.

    distfs = {'euclidean', 'squaredeuclidean', 'seuclidean'};
    if any(1 == size(data)), data = data(:); end
    if size(data, 2) > 1
        distfs = [distfs, {'cosine', 'correlation', 'spearman', 'jaccard', ...
                'mahalanobis'}];
    end
    
    if size(data, 1) < 2
        error('rsa_algorithms.create_rsm: Not enough observations');
    end
    
    if iscell(data)
        if ~ischar(data{1}) 
            error(['Invalid input to rsa_algorithms.create_rsm: ', ...
                'string cells expected']); 
        end
        rsm = eye(size(data,1)) + squareform(pdist(data, @featuresim));
        if nargin>1 
            disp('Warning: distance function directive ignored'); end
    elseif isinteger(data)
        rsm = eye(size(data,1)) + squareform(pdist(data, @categorysim));
        if nargin>1
            disp('Warning: distance function directive ignored'); end
    elseif isfloat(data)
        if (ischar(distf) && ~ismember({distf}, distfs)) && ...
                ~isa(distf, 'function_handle')
            error(['Invalid input to rsa_algorithms.create_rsm: ', ...
                'invalid distance function']);
        end
        if ischar(distf) && strcmp(distf, 'correlation')
            rsm = 1 - squareform(pdist(data, distf));
        elseif ischar(distf)
            %since a correlation will be used, the conversion to similarity
            %will be ok irrespective of scale and location.
            rsm = squareform(pdist(data, distf));
            rsm = 1 - rsm ./ max(rsm(:));
        else   %function handle
            rsm = squareform(pdist(data, distf));
        end
    else
        error(['rsa_algorithms.create_rsm: Do not know how to use ',...
            class(data), ' type in data']);
    end
    
    if any(isnan(rsm(:)))
        error('rsa_algorithms.create_rsm: invalid data (NaNs generated)');
    end
    
    function value = featuresim(X, XX)
        value = zeros(1, size(XX, 1));
        for i = 1 : size(XX, 1)
            value(i) = sum(cellfun(@strcmp, X, XX(i,:)))./numel(X);
        end
    end

    function value = categorysim(X, XX)
        value = zeros(1, size(XX, 1));
        for i = 1 : size(XX, 1)
            value(i) = sum(X == XX(i,:)) ./ numel(X);
        end
    end    
end

