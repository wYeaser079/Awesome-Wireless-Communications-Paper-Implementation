%RUN_ALL Batch-run main_figure4/5/6 with reduced Monte-Carlo budget.
%   Also patches functionComputeSE_uplink/downlink to print per-realization
%   progress, and patches the three main scripts to use a small
%   nbrOfSetups/nbrOfRealizations so the run fits in a laptop session.
%   Originals are restored on cleanup.

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
out_dir = fullfile(here, 'outputs');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

main_scripts = {'main_figure4', 'main_figure5', 'main_figure6'};
save_names = {
    {'fig4a', 'fig4b'};
    {'fig5a', 'fig5b'};
    {'fig6a', 'fig6b'};
    };
helper_scripts = {'functionComputeSE_uplink', 'functionComputeSE_downlink'};

% --- Back up every file we touch, then patch ---
all_files = [main_scripts, helper_scripts];
backups = cell(length(all_files), 2);
for i = 1:length(all_files)
    p_path = fullfile(here, [all_files{i} '.m']);
    orig = fileread(p_path);
    backups{i,1} = p_path;
    backups{i,2} = orig;
end

% Patch main scripts: reduce Monte-Carlo budget
for i = 1:length(main_scripts)
    p_path = backups{i,1};
    orig = backups{i,2};
    patched = regexprep(orig, 'nbrOfSetups\s*=\s*\d+\s*;', ...
                              'nbrOfSetups = 2;');
    patched = regexprep(patched, 'nbrOfRealizations\s*=\s*\d+\s*;', ...
                                 'nbrOfRealizations = 20;');
    fid = fopen(p_path, 'w'); fwrite(fid, patched); fclose(fid);
end

% Patch helper files: insert per-realization progress print after
% each `for nreal = 1:nbrOfRealizations` line.
for i = length(main_scripts)+1 : length(all_files)
    p_path = backups{i,1};
    orig = backups{i,2};
    progress_line = ['    if mod(nreal, 2) == 0 || nreal == 1, ' ...
                     'fprintf(''     nreal %d/%d\\n'', nreal, nbrOfRealizations); end'];
    patched = regexprep(orig, ...
        '(for nreal = 1:nbrOfRealizations)', ...
        ['$1' char(10) progress_line]);
    fid = fopen(p_path, 'w'); fwrite(fid, patched); fclose(fid);
end

cleanupObj = onCleanup(@() restore_all(backups));
rehash;

master_log = fullfile(out_dir, 'run_all_master.log');
mfid = fopen(master_log, 'w');
fprintf(mfid, 'RUN_ALL master log\nStarted: %s\n\n', datestr(now)); %#ok<*DATST,*TNOW1>

for j = 1:length(main_scripts)
    script = main_scripts{j};
    names  = save_names{j};
    log_path = fullfile(out_dir, [script '.log']);

    fprintf('\n============================================================\n');
    fprintf('RUNNING: %s\n', script);
    fprintf('============================================================\n');
    fprintf(mfid, '\n=== %s ===\n', script);

    t0 = tic;
    try
        diary(log_path); diary on;
        close all;
        run_script(script);
        diary off;

        figs = findobj('Type', 'figure');
        [~, ord] = sort([figs.Number]);
        figs = figs(ord);
        for i = 1:numel(figs)
            f = figs(i);
            if i <= numel(names)
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

fprintf(mfid, '\nFinished: %s\n', datestr(now));
fclose(mfid);
fprintf('\nAll done. Outputs in: %s\n', out_dir);


function run_script(script_name)
    run(script_name);
end


function restore_all(backups)
    for i = 1:size(backups, 1)
        fid = fopen(backups{i,1}, 'w');
        fwrite(fid, backups{i,2});
        fclose(fid);
    end
    fprintf('main and helper scripts restored.\n');
end
