function rsm = rsa_create_rsm(data, name, distf)
% RSA_CREATE_RSM creates a representational similarity structure for a
% matrix of data, where each row is an observation and each column the
% property values. Name gives the name of the object class characterized by
% the similarity.
%
% Arguments
% data can contain the following types:
%   a cell matrix of strings: each column is the value of a qualitative
%       property
%   a matrix of integers: (convert numbers with int16(), uint16(), etc. to
%       form a matrix of integers) each number value is treated as a
%       qualitative factor level
%   a matrix of double values: computes similarity according to distf,
%       which is one of ''cosine'', ''correlation'', ''spearman'',
%       ''jaccard''. If distf is not provided or empty, defaults to
%       ''correlation'' or ''euclidean'' if just one column. See Matlab's 
%       pdist function for an explanation of these distance functions. If 
%       data is not double, distf is ignored.
% name: the name of the map
% distf:the distance function used (see explanation for data argument).
%
% Output
% The rsm structure contains the following fields:
%   name: the name of the map
%   type: 'similarity'
%   RSM:  a matrix with the similarity values. Only the upper triangular
%         portion is ever used.
%   data: the data used to compute the RSM (not in use at present)
%   distf:the name of the distance function used (not in use at present)
%
% 2020-21 Roberto Viviani - December 2020
% Institute of Psychology, University of Innsbruck
% Last modified: August 1st 2021


    if nargin < 2
        error('rsa_create_rsm: Not enough arguments: name required.');
    end
    if nargin < 1
        error('rsa_create_rsm: Not enough arguments: features and name required');
    end
    
    %third argument will be used only if data is composed of floats.
    if nargin < 3 || isempty(distf)
        if size(data, 2) == 1
            distf = 'euclidean';
        else
            distf = 'correlation'; 
        end
    end

    rsaio = rsa_rsmio();
    rsaalgo = rsa_algorithms();
    if iscell(data) || isinteger(data)
        rsm = rsaio.create_rsm(name, rsaalgo.create_rsm(data));
        distf = 'categories';
        if nargin>2, disp('Warning: distance function directive ignored'); end
    elseif isfloat(data)
        rsm = rsaio.create_rsm(name, rsaalgo.create_rsm(data, distf));
    else
        error(['rsa_create_rsm: Do not know how to use type ', class(data)]);
    end
    rsm.data = data;
    rsm.distf = distf;
end


