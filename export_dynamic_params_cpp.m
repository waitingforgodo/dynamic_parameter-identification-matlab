% Export identified pi_b, pi_fr to cpp/include/dynamic_params.h
% Run after estimate_dynamic_params in main.m:
%   export_dynamic_params_cpp(sol)

function export_dynamic_params_cpp(sol, outFile)
if nargin < 2 || isempty(outFile)
    outFile = fullfile('cpp_le', 'include', 'dynamic_params.h');
end
writeVectorH(outFile, 'PI_B', sol.pi_b, 'PI_FR', sol.pi_fr);
fprintf('Exported dynamic params to %s\n', outFile);
end

function writeVectorH(path, nameB, pi_b, nameFr, pi_fr)
pi_b = pi_b(:);
pi_fr = pi_fr(:);
fid = fopen(path, 'w');
fprintf(fid, '#pragma once\n\n');
fprintf(fid, '#include "robot_dyn_config.h"\n\n');
fprintf(fid, 'static const double %s[N_BASE_PARAMS] = {\n', nameB);
for k = 1:numel(pi_b)
    fprintf(fid, '    %.17e%s\n', pi_b(k), ternary(k < numel(pi_b), ',', ''));
end
fprintf(fid, '};\n\n');
fprintf(fid, 'static const double %s[N_FRICTION_PARAMS] = {\n', nameFr);
for k = 1:numel(pi_fr)
    fprintf(fid, '    %.17e%s\n', pi_fr(k), ternary(k < numel(pi_fr), ',', ''));
end
fprintf(fid, '};\n');
fclose(fid);
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
