function start_parallel_pool(worker_count)
%START_PARALLEL_POOL Start a local parallel pool when none is active.

if nargin < 1 || isempty(worker_count)
    worker_count = feature('numcores');
end

pool = gcp('nocreate');
if isempty(pool)
    parpool('local', worker_count);
end
end

