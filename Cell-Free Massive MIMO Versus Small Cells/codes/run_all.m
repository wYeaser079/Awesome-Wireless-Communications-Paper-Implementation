%RUN_ALL Batch-run every main_*.m script, save figures to PNG, capture text.
%   Runs each script in its own function workspace, saves all figures
%   produced (numbered fig1.png, fig2.png, ...) into the outputs/ folder,
%   and writes a per-script log file.
%
%   To keep total runtime reasonable on a laptop, a reduced num_setups is
%   injected by temporarily overriding params.m on the MATLAB path.

clear; clc; close all;

out_dir = fullfile(fileparts(mfilename('fullpath')), 'outputs');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% Back up and patch params.m to reduce Monte Carlo load
params_path = fullfile(fileparts(mfilename('fullpath')), 'params.m');
backup_path = fullfile(fileparts(mfilename('fullpath')), 'params_backup.m.bak');
orig_params = fileread(params_path);
fid = fopen(backup_path, 'w'); fwrite(fid, orig_params); fclose(fid);

patched = regexprep(orig_params, ...
    'p\.num_setups\s*=\s*\d+\s*;', 'p.num_setups = 15;');
patched = regexprep(patched, ...
    'p\.bisection_maxiter\s*=\s*\d+\s*;', 'p.bisection_maxiter = 15;');
fid = fopen(params_path, 'w'); fwrite(fid, patched); fclose(fid);
cleanupObj = onCleanup(@() restore_params(params_path, orig_params));
rehash;

% Scripts and the figure base-names to save under
jobs = {
    'main_fig2',          {'fig2'};
    'main_fig9',          {'fig9'};
    'main_fig10',         {'fig10'};
    'main_fig11',         {'fig11'};
    'main_fig3_fig4',     {'fig3','fig4'};
    'main_fig5_fig6',     {'fig5','fig6'};
    'main_fig7_fig8',     {'fig7','fig8'};
    'main_table2_table3', {}
    };

master_log_path = fullfile(out_dir, 'run_all_master.log');
mfid = fopen(master_log_path, 'w');
fprintf(mfid, 'RUN_ALL master log\nStarted: %s\n\n', datestr(now)); %#ok<DATST,TNOW1>

for j = 1:size(jobs,1)
    script = jobs{j,1};
    names  = jobs{j,2};
    log_path = fullfile(out_dir, [script '.log']);

    fprintf('\n============================================================\n');
    fprintf('RUNNING: %s\n', script);
    fprintf('============================================================\n');
    fprintf(mfid, '\n=== %s ===\n', script);

    t0 = tic;
    try
        diary(log_path); diary on;
        close all;
        run_script(script);   % run in helper so script's "clear" is scoped
        diary off;

        % Save all open figures
        figs = findobj('Type', 'figure');
        figs = flipud(figs);  % preserve creation order
        for i = 1:numel(figs)
            f = figs(i);
            if i <= numel(names) && ~isempty(names)
                base = names{i};
            else
                base = sprintf('%s_fig%d', script, i);
            end
            png_path = fullfile(out_dir, [base '.png']);
            try
                exportgraphics(f, png_path, 'Resolution', 150);
            catch
                saveas(f, png_path);
            end
            fprintf('  Saved %s\n', png_path);
        end
        close all;
        elapsed = toc(t0);
        fprintf(mfid, '  SUCCESS (%.1f s)\n', elapsed);
    catch ME
        diary off;
        elapsed = toc(t0);
        fprintf(2, 'ERROR in %s: %s\n', script, ME.message);
        fprintf(mfid, '  FAILED (%.1f s): %s\n', elapsed, ME.message);
        fprintf(mfid, '  %s\n', getReport(ME, 'basic'));
    end
end

fprintf(mfid, '\nFinished: %s\n', datestr(now)); %#ok<DATST,TNOW1>
fclose(mfid);
fprintf('\nAll done. Outputs in: %s\n', out_dir);


function run_script(script_name)
    % Run a script that starts with 'clear; clc; close all;'.
    % Calling run() inside this helper makes the script share this
    % helper's local workspace, so the script's 'clear' is scoped.
    run(script_name);
end


function restore_params(params_path, orig)
    fid = fopen(params_path, 'w');
    fwrite(fid, orig);
    fclose(fid);
    fprintf('params.m restored.\n');
end
