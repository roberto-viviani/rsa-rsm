function args = ask_args(varargin)
% ASK_ARGS This is a helper function to setup function calls
% of functions with prefix exp_ so that the function does not
% start, but returns an argument structure, or, when called with
% an argument structure, asks user to modify specific parameters,
% without executing the function
%
% Usage examples (using exp_t)
%
%   args = exp_t(ask_args);             %ask arguments for exp_t
%   args = exp_t(ask_args(args));       %ask arguments missing in args
%           %asks new pars for ClusterTest and OutputDir
%   args = exp_t(ask_args(args, 'ClusterTest', 'OutputDir'));
%
% 2009 Roberto Viviani 
% Abteilung Psychiatrie III, Universitätsklinikum Ulm
% Last modified: 20 august 2009

    switch nargin
    case 0
        args = struct('Args', 1);
    case 1
        args = varargin{1};
        args.Args = 1;
    otherwise
        args = varargin{1};
        for a = 2 : nargin
            if isfield(args, varargin{a})
                args = rmfield(args, varargin{a});
            end
        end
        args.Args = 1;
    end
         
