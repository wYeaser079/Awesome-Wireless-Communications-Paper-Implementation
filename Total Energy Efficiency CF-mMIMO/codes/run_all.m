%RUN_ALL Batch-run main_Fig1..6, save figures and per-script logs.
%   The defaults already use small (M,K,NMC) grids, so no patching is
%   needed. We just dispatch each script, capture its figure, save to
%   outputs/figN.png, and write a per-script diary.

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
out_dir = fullfile(here, 'outputs');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

jobs = {
    'main_Fig1_convergence',      {'fig1'};
    'main_Fig2_power_control',    {'fig2'};
    'main_Fig3_AP_selection_Pbt', {'fig3'};
    'main_Fig4_num_APs_per_user', {'fig4'};
    'main_Fig5_effect_of_N',      {'fig5'};
    'main_Fig6_CF_vs_colocated',  {'fig6'};
    };

master_log = fullfile(out_dir, 'run_all_master.log');
mfid = fopen(master_log, 'w');
fprintf(mfid, 'RUN_ALL master log\nStarted: %s\n\n', datestr(now)); %#ok<*DATST,*TNOW1>

for j = 1:size(jobs, 1)
    script = jobs{j, 1};
    names  = jobs{j, 2};
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
