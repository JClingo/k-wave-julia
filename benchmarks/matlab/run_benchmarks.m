%% KWave MATLAB Performance Benchmark Suite
%  ==========================================
%  Benchmarks full simulations (1D / 2D / 3D), FFT throughput, and GPU
%  performance using the k-Wave MATLAB toolbox.
%
%  Prerequisites:
%    - k-Wave toolbox on the MATLAB path  (http://www.k-wave.org) or cloned from GitHub or installed via Add-Ons
%    - MATLAB R2021a or later
%    - Parallel Computing Toolbox  (for GPU section)
%
%  Usage:
%    run_benchmarks          % CPU only
%    run_benchmarks('gpu')   % CPU + GPU
%
%  Results written to:
%    benchmarks/results/matlab/full_simulations.csv
%    benchmarks/results/matlab/fft_components.csv
%    benchmarks/results/matlab/gpu_simulations.csv  (if GPU requested)

function run_benchmarks(varargin)

    do_gpu = ~isempty(varargin) && strcmpi(varargin{1}, 'gpu');

    %% Config
    N_WARMUP  = 1;
    N_SAMPLES = 5;
    DX        = 1e-4;   % 0.1 mm
    C0        = 1500;   % m/s
    RHO0      = 1000;   % kg/m³
    RNG_SEED  = 42;

    results_dir = fullfile(fileparts(mfilename('fullpath')), ...
                           '..', 'results', 'matlab');
    if ~exist(results_dir, 'dir'); mkdir(results_dir); end

    fprintf('=================================================================\n');
    fprintf('k-Wave MATLAB Benchmark Suite\n');
    fprintf('MATLAB:    %s\n', version);
    try
        kw_ver = kWaveVersion;
        fprintf('k-Wave:    %s\n', kw_ver);
    catch
        fprintf('k-Wave:    version unavailable\n');
    end
    fprintf('Date:      %s\n', datestr(now));
    fprintf('N_WARMUP=%d  N_SAMPLES=%d\n', N_WARMUP, N_SAMPLES);
    fprintf('=================================================================\n');

    %% ── Check k-Wave is on path ─────────────────────────────────────────
    if ~exist('kspaceFirstOrder1D', 'file')
        error(['k-Wave toolbox not found on MATLAB path.\n' ...
               'Add k-Wave root to path: addpath(''/path/to/k-wave'')']);
    end

    %% ── Scenario definitions ────────────────────────────────────────────
    %  {label, dim, Nx, Ny, Nz, t_end, pml_size, medium_type}
    scenarios = {
        % 1D
        '1d_small',       1,  256,   1,  1, 4e-6, 20, 'homogeneous';
        '1d_medium',      1, 1024,   1,  1, 4e-6, 20, 'homogeneous';
        '1d_large',       1, 4096,   1,  1, 4e-6, 20, 'homogeneous';
        '1d_absorbing',   1, 1024,   1,  1, 4e-6, 20, 'absorbing';
        '1d_hetero',      1, 1024,   1,  1, 4e-6, 20, 'heterogeneous';
        % 2D
        '2d_small',       2,   64,  64,  1, 4e-6, 20, 'homogeneous';
        '2d_medium',      2,  256, 256,  1, 4e-6, 20, 'homogeneous';
        '2d_large',       2,  512, 512,  1, 2e-6, 20, 'homogeneous';
        '2d_absorbing',   2,  256, 256,  1, 4e-6, 20, 'absorbing';
        '2d_hetero',      2,  256, 256,  1, 4e-6, 20, 'heterogeneous';
        % 3D
        '3d_small',       3,   32,  32, 32, 2e-6,  8, 'homogeneous';
        '3d_medium',      3,   64,  64, 64, 2e-6, 12, 'homogeneous';
        '3d_large',       3,  128, 128,128, 1e-6, 20, 'homogeneous';
        '3d_absorbing',   3,   64,  64, 64, 2e-6, 12, 'absorbing';
        '3d_hetero',      3,   64,  64, 64, 2e-6, 12, 'heterogeneous';
    };

    %% ── CPU Float64 benchmarks ──────────────────────────────────────────
    fprintf('\n─── Full-simulation benchmarks  backend=cpu_f64 ───\n');
    sim_rows_f64 = run_sim_scenarios(scenarios, N_WARMUP, N_SAMPLES, ...
                                     DX, C0, RHO0, RNG_SEED, 'cpu_f64', 'double');

    %% ── CPU Single-precision benchmarks ────────────────────────────────
    fprintf('\n─── Full-simulation benchmarks  backend=cpu_f32 ───\n');
    % homogeneous only (matches Julia cpu_f32 set)
    homog_mask = strcmp(scenarios(:,8), 'homogeneous');
    sim_rows_f32 = run_sim_scenarios(scenarios(homog_mask,:), N_WARMUP, N_SAMPLES, ...
                                     DX, C0, RHO0, RNG_SEED, 'cpu_f32', 'single');

    % Combine and save full-simulation results
    sim_rows = [sim_rows_f64, sim_rows_f32];
    T = struct2table(vertcat(sim_rows{:}));
    writetable(T, fullfile(results_dir, 'full_simulations.csv'));
    fprintf('\nSaved: full_simulations.csv\n');

    %% ── FFT component benchmarks ────────────────────────────────────────
    fprintf('\n─── FFT component micro-benchmarks ───\n');
    fft_rows = run_fft_benchmarks();
    T_fft = struct2table(vertcat(fft_rows{:}));
    writetable(T_fft, fullfile(results_dir, 'fft_components.csv'));
    fprintf('Saved: fft_components.csv\n');

    %% ── GPU benchmarks (optional) ───────────────────────────────────────
    if do_gpu
        gpu_scenarios = {
            '2d_small',   2,  64,  64,  1, 4e-6, 20, 'homogeneous';
            '2d_medium',  2, 256, 256,  1, 4e-6, 20, 'homogeneous';
            '2d_large',   2, 512, 512,  1, 2e-6, 20, 'homogeneous';
            '3d_small',   3,  32,  32, 32, 2e-6,  8, 'homogeneous';
            '3d_medium',  3,  64,  64, 64, 2e-6, 12, 'homogeneous';
            '3d_large',   3, 128, 128,128, 1e-6, 20, 'homogeneous';
        };

        fprintf('\n─── GPU benchmarks  backend=gpu_f32 ───\n');
        % Check GPU availability
        try
            gpu_info = gpuDevice;
            fprintf('GPU: %s\n', gpu_info.Name);
        catch
            fprintf('No GPU available — skipping GPU benchmarks.\n');
            fprintf('Done.  %s\n', datestr(now));
            return;
        end

        gpu_rows = run_sim_scenarios(gpu_scenarios, N_WARMUP, N_SAMPLES, ...
                                     DX, C0, RHO0, RNG_SEED, 'gpu_f32', 'gpuArray-single');
        T_gpu = struct2table(vertcat(gpu_rows{:}));
        writetable(T_gpu, fullfile(results_dir, 'gpu_simulations.csv'));
        fprintf('Saved: gpu_simulations.csv\n');
    end

    fprintf('\nDone.  %s\n', datestr(now));
end

%% ──────────────────────────────────────────────────────────────────────────
%  run_sim_scenarios — run a set of scenarios and return result rows
%% ──────────────────────────────────────────────────────────────────────────
function rows = run_sim_scenarios(scenarios, N_WARMUP, N_SAMPLES, ...
                                  DX, C0, RHO0, RNG_SEED, backend, data_cast)
    rows = {};
    rng(RNG_SEED);   % reproducible heterogeneous maps

    for s_idx = 1:size(scenarios, 1)
        label      = scenarios{s_idx, 1};
        dim        = scenarios{s_idx, 2};
        Nx         = scenarios{s_idx, 3};
        Ny         = scenarios{s_idx, 4};
        Nz         = scenarios{s_idx, 5};
        t_end      = scenarios{s_idx, 6};
        pml_size   = scenarios{s_idx, 7};
        mtype      = scenarios{s_idx, 8};

        % ── Build kgrid ──────────────────────────────────────────────────
        switch dim
            case 1
                kgrid = kWaveGrid(Nx, DX);
            case 2
                kgrid = kWaveGrid(Nx, DX, Ny, DX);
            case 3
                kgrid = kWaveGrid(Nx, DX, Ny, DX, Nz, DX);
        end
        kgrid.makeTime(C0, 0.3, t_end);
        Nt = kgrid.Nt;

        % ── Build medium ─────────────────────────────────────────────────
        switch mtype
            case 'homogeneous'
                medium.sound_speed = C0;
                medium.density     = RHO0;
            case 'absorbing'
                medium.sound_speed = C0;
                medium.density     = RHO0;
                medium.alpha_coeff = 0.75;
                medium.alpha_power = 1.5;
            case 'heterogeneous'
                switch dim
                    case 1
                        c_map = C0   * (1 + 0.10 * randn(Nx, 1));
                        rho_m = RHO0 * ones(Nx, 1);
                    case 2
                        c_map = C0   * (1 + 0.10 * randn(Nx, Ny));
                        rho_m = RHO0 * ones(Nx, Ny);
                    case 3
                        c_map = C0   * (1 + 0.10 * randn(Nx, Ny, Nz));
                        rho_m = RHO0 * ones(Nx, Ny, Nz);
                end
                medium.sound_speed = c_map;
                medium.density     = rho_m;
        end

        % ── Build source ─────────────────────────────────────────────────
        switch dim
            case 1
                p0 = zeros(Nx, 1);
                p0(Nx/2-5:Nx/2+5) = 1;
                source.p0 = p0;
            case 2
                source.p0 = makeDisc(Nx, Ny, Nx/2, Ny/2, 5);
            case 3
                source.p0 = makeBall(Nx, Ny, Nz, Nx/2, Ny/2, Nz/2, 3);
        end

        % ── Build sensor ─────────────────────────────────────────────────
        switch dim
            case 1
                sensor.mask = zeros(Nx, 1);
                sensor.mask(round(Nx/4))   = 1;
                sensor.mask(round(3*Nx/4)) = 1;
            case 2
                sensor.mask = zeros(Nx, Ny);
                sensor.mask(round(Nx/4),   round(Ny/2)) = 1;
                sensor.mask(round(3*Nx/4), round(Ny/2)) = 1;
            case 3
                sensor.mask = zeros(Nx, Ny, Nz);
                sensor.mask(round(Nx/4),   round(Ny/2), round(Nz/2)) = 1;
                sensor.mask(round(3*Nx/4), round(Ny/2), round(Nz/2)) = 1;
        end

        % ── Extra options ────────────────────────────────────────────────
        extra_opts = {'Smooth', false, 'PlotSim', false, ...
                      'PMLSize', pml_size};
        if ~strcmp(data_cast, 'double')
            extra_opts = [extra_opts, 'DataCast', data_cast]; %#ok<AGROW>
        end

        % ── Solver function handle ────────────────────────────────────────
        switch dim
            case 1; solver = @() kspaceFirstOrder1D(kgrid, medium, source, sensor, extra_opts{:});
            case 2; solver = @() kspaceFirstOrder2D(kgrid, medium, source, sensor, extra_opts{:});
            case 3; solver = @() kspaceFirstOrder3D(kgrid, medium, source, sensor, extra_opts{:});
        end

        fprintf('  %-22s  Nx=%-5d Nt=%-5d  ', label, Nx, Nt);

        try
            r = bench_fn(solver, N_WARMUP, N_SAMPLES);

            fprintf('first=%6.3fs  median=%6.3fs\n', r.first_run_s, r.median_s);

            row.timestamp    = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
            row.language     = 'matlab';
            row.backend      = backend;
            row.scenario     = label;
            row.dim          = dim;
            row.Nx           = Nx;
            row.Ny           = Ny;
            row.Nz           = Nz;
            row.Nt           = Nt;
            row.medium_type  = mtype;
            row.first_run_s  = r.first_run_s;
            row.min_s        = r.min_s;
            row.median_s     = r.median_s;
            row.mean_s       = r.mean_s;
            row.max_s        = r.max_s;
            row.std_s        = r.std_s;
            row.median_alloc = 0;    % not tracked in MATLAB
            row.median_gc_s  = 0;

            rows{end+1} = row; %#ok<AGROW>
        catch ME
            fprintf(' ERROR: %s\n', ME.message);
        end
    end
end

%% ──────────────────────────────────────────────────────────────────────────
%  run_fft_benchmarks — standalone FFT throughput
%% ──────────────────────────────────────────────────────────────────────────
function rows = run_fft_benchmarks()
    rows = {};
    N_FFT_SAMPLES = 500;

    % 1D
    for N = [128, 256, 512, 1024, 2048, 4096, 8192, 16384]
        x = rand(N, 1) + 1i*rand(N, 1);
        % warmup
        for j = 1:5; fft(x); end
        times = zeros(N_FFT_SAMPLES, 1);
        for j = 1:N_FFT_SAMPLES
            t = tic; fft(x); times(j) = toc(t);
        end
        times_ns  = times * 1e9;
        flops     = 5 * N * log2(N);
        med_ns    = median(times_ns);
        fprintf('  FFT 1D  N=%-6d  median=%7.2f µs  %.2f GFLOP/s\n', ...
                N, med_ns/1e3, flops/med_ns);

        row.language = 'matlab'; row.backend = 'cpu_f64';
        row.label = 'fft_1d'; row.dim = 1;
        row.N = N; row.Nx = N; row.Ny = 1; row.Nz = 1;
        row.median_ns = med_ns; row.min_ns = min(times_ns);
        row.throughput_gflops = flops / med_ns;
        rows{end+1} = row; %#ok<AGROW>
    end

    % 2D
    sizes_2d = [64, 128, 256, 512, 1024];
    for k = 1:length(sizes_2d)
        Nx = sizes_2d(k); Ny = Nx;
        x = rand(Nx, Ny) + 1i*rand(Nx, Ny);
        for j = 1:3; fft2(x); end
        times = zeros(100, 1);
        for j = 1:100
            t = tic; fft2(x); times(j) = toc(t);
        end
        times_ns  = times * 1e9;
        N_total   = Nx * Ny;
        flops     = 5 * N_total * log2(N_total);
        med_ns    = median(times_ns);
        fprintf('  FFT 2D  %4dx%-4d  median=%7.2f ms  %.2f GFLOP/s\n', ...
                Nx, Ny, med_ns/1e6, flops/med_ns);

        row.language = 'matlab'; row.backend = 'cpu_f64';
        row.label = 'fft_2d'; row.dim = 2;
        row.N = N_total; row.Nx = Nx; row.Ny = Ny; row.Nz = 1;
        row.median_ns = med_ns; row.min_ns = min(times_ns);
        row.throughput_gflops = flops / med_ns;
        rows{end+1} = row; %#ok<AGROW>
    end

    % 3D
    sizes_3d = [32, 64, 128, 256];
    for k = 1:length(sizes_3d)
        Nx = sizes_3d(k); Ny = Nx; Nz = Nx;
        x = rand(Nx, Ny, Nz) + 1i*rand(Nx, Ny, Nz);
        for j = 1:3; fftn(x); end
        times = zeros(50, 1);
        for j = 1:50
            t = tic; fftn(x); times(j) = toc(t);
        end
        times_ns  = times * 1e9;
        N_total   = Nx * Ny * Nz;
        flops     = 5 * N_total * log2(N_total);
        med_ns    = median(times_ns);
        fprintf('  FFT 3D  %dx%dx%-3d  median=%7.2f ms  %.2f GFLOP/s\n', ...
                Nx, Ny, Nz, med_ns/1e6, flops/med_ns);

        row.language = 'matlab'; row.backend = 'cpu_f64';
        row.label = 'fft_3d'; row.dim = 3;
        row.N = N_total; row.Nx = Nx; row.Ny = Ny; row.Nz = Nz;
        row.median_ns = med_ns; row.min_ns = min(times_ns);
        row.throughput_gflops = flops / med_ns;
        rows{end+1} = row; %#ok<AGROW>
    end
end

%% ──────────────────────────────────────────────────────────────────────────
%  bench_fn — warmup + N_SAMPLES timed calls of a function handle
%% ──────────────────────────────────────────────────────────────────────────
function r = bench_fn(f, N_WARMUP, N_SAMPLES)
    % First call — includes any lazy initialisation overhead
    t0 = tic; f(); r.first_run_s = toc(t0);

    % Additional warmup
    for i = 2:N_WARMUP; f(); end

    % Timed samples
    times = zeros(N_SAMPLES, 1);
    for i = 1:N_SAMPLES
        t0 = tic; f(); times(i) = toc(t0);
    end

    r.min_s    = min(times);
    r.median_s = median(times);
    r.mean_s   = mean(times);
    r.max_s    = max(times);
    r.std_s    = std(times);
end
