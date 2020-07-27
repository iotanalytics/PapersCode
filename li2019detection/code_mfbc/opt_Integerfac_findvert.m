% Function to set-up options for 'Integerfac_findvert.m'
%
% Inputs
%
% 'affine'  --- Are sum-to-one constraints imposed on A (default: true) 
%
% 'nonnegative' --- Are non-negativity constraints imposed on A (default: true)
%
%
% 'aggr'    --- One of 'no' (default), 'bestsingle', 'backward'
%                 
%
% 'replace' --- number between 0 and 1. Gives the fraction of 
%               coordinates to be re-used after each coordinate
%               selection. If 'replace = 1' (default) all
%               coordinates are discarded from subsequent row
%               selections. 
%
% 'nsamples' --- number of different coordinate choices
%                          to be considered. 
%                          '0': use default setting
%                          integer: number of samples to be considered
%                          negative integer: extra samples (chosen
%                          randomly) to be considered.
%
%
% 'naggsets' --- the naggsets top best single fits are considered for
%                   aggregation (default: 5).
%
% 'chunksize' --- chunksize for splitting the right hand sides
%                  occuring in the solutions of the linear systems
%                  (default: 10) 
% 
%
% 'verbose' --- status is information is displayed/notdisplayed (true[default]/false).
%
%
% Output:
%
%
% A struct to be passed on as argument 'opt' to 'Integerfac_findvert.m'
             

function options = opt_Integerfac_findvert(varargin)

p = inputParser;


valid_affine = @(x) validateattributes(x, {'logical'}, {'scalar'});
p.addParamValue('affine', true, valid_affine);

valid_nonnegative = @(x) validateattributes(x, {'logical'}, {'scalar'});
p.addParamValue('nonnegative', true, valid_nonnegative);

valid_aggr = @(x) ismember(x, {'no', 'bestsingle', 'backward'});
p.addParamValue('aggr', 'no', valid_aggr);

valid_replace = @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', '<=', 1});
p.addParamValue('replace', 1, valid_replace);

valid_nsamples = @(x) validateattributes(x,{'numeric'},{'scalar','integer'});
p.addParamValue('nsamples', 0, valid_nsamples);

valid_naggsets = @(x) validateattributes(x,{'numeric'},{'scalar','integer'});
p.addParamValue('naggsets', 5, valid_naggsets);

valid_chunksize = @(x) validateattributes(x,{'numeric'},{'scalar','integer'});
p.addParamValue('chunksize', 10, valid_chunksize);

valid_verbose = @(x) validateattributes(x, {'logical'}, {'scalar'});
p.addParamValue('verbose', true, valid_verbose);

p.parse(varargin{:});
options = p.Results;

