# Shared shell functions for the environment tests. Source this file from an sbatch
# script after `module load Miniforge3/...`.
#
# Everything (environments, package caches, HOME, TMPDIR) lives under WORK, a folder on
# the node-local disk /tmp that disappears with the job.

# isolate WORK
#   Points HOME, TMPDIR and every conda, pip, Jupyter, IPython, matplotlib and numba
#   cache to folders under WORK, and removes interpreter settings of the login shell.
#   A fresh HOME also hides ~/.condarc, ~/.config/pip, ~/.jupyter, ~/.ipython and the
#   user site ~/.local/lib/pythonX.Y.
isolate() {
    local work="$1"
    mkdir -p "$work"/{home,tmp,pkgs,pipcache,cache,envs,mamba}
    export HOME="$work/home"
    export TMPDIR="$work/tmp"
    export XDG_CACHE_HOME="$work/cache"
    export XDG_CONFIG_HOME="$work/home/.config"
    export XDG_DATA_HOME="$work/home/.local/share"
    export CONDA_PKGS_DIRS="$work/pkgs"
    export CONDA_ENVS_DIRS="$work/envs"
    export MAMBA_ROOT_PREFIX="$work/mamba"
    export CONDA_REGISTER_ENVS=false
    export CONDA_CHANNEL_PRIORITY=strict
    export PIP_CACHE_DIR="$work/pipcache"
    export PIP_CONFIG_FILE=/dev/null
    export PIP_DISABLE_PIP_VERSION_CHECK=1
    export JUPYTER_CONFIG_DIR="$work/home/.jupyter"
    export JUPYTER_DATA_DIR="$work/home/.local/share/jupyter"
    export JUPYTER_RUNTIME_DIR="$work/home/.local/share/jupyter/runtime"
    export IPYTHONDIR="$work/home/.ipython"
    export MPLCONFIGDIR="$work/home/.config/matplotlib"
    export NUMBA_CACHE_DIR="$work/cache/numba"
    unset PYTHONPATH PYTHONHOME PYTHONSTARTUP PYTHONUSERBASE
    unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL CONDA_PROMPT_MODIFIER CONDA_EXE CONDA_PYTHON_EXE _CE_M _CE_CONDA
    unset R_LIBS R_LIBS_USER
    export PYTHONNOUSERSITE=1
}

# timed LABEL LOGFILE cmd...
#   Runs a command under GNU time and appends "LABEL elapsed=... maxrss_kb=..." to
#   LOGFILE. The command's own output goes to stdout/stderr as usual.
timed() {
    local label="$1" logfile="$2"; shift 2
    local tfile="$TMPDIR/time_${label}.txt"
    local rc=0
    /usr/bin/time -v -o "$tfile" "$@" || rc=$?
    local el rss
    el=$(grep "Elapsed (wall clock)" "$tfile" | awk '{print $NF}')
    rss=$(grep "Maximum resident set size" "$tfile" | awk '{print $NF}')
    echo "$label exit=$rc elapsed=$el maxrss_kb=$rss" | tee -a "$logfile"
    return $rc
}

# record_env PREFIX NAME STAGING LOGDIR
#   Writes the four records of an environment:
#     STAGING/NAME.lock.txt          conda list --explicit --md5
#     STAGING/NAME.pip-freeze.txt    pip list --format=freeze (name==version of every
#                                    Python distribution, conda- or pip-installed)
#     LOGDIR/NAME.conda-env-export.yml   conda env export --no-builds
#     LOGDIR/NAME.conda-list.txt         conda list (shows the channel of each package)
#     LOGDIR/NAME.pip-freeze-raw.txt     pip freeze as printed
#     LOGDIR/NAME.python-version.txt     sys.version
record_env() {
    local prefix="$1" name="$2" staging="$3" logdir="$4"
    local conda="$EBROOTMINIFORGE3/bin/conda"
    "$conda" list -p "$prefix" --explicit --md5 > "$staging/$name.lock.txt"
    "$conda" env export -p "$prefix" --no-builds | sed '/^prefix:/d' | sed "s/^name: .*/name: $name/" > "$logdir/$name.conda-env-export.yml"
    "$conda" list -p "$prefix" > "$logdir/$name.conda-list.txt"
    "$prefix/bin/python" -m pip list --format=freeze > "$staging/$name.pip-freeze.txt"
    "$prefix/bin/python" -m pip freeze > "$logdir/$name.pip-freeze-raw.txt"
    "$prefix/bin/python" -c "import sys; print(sys.version)" > "$logdir/$name.python-version.txt"
    echo "records written for $name: $(wc -l < "$staging/$name.lock.txt") lock lines, $(wc -l < "$staging/$name.pip-freeze.txt") pip lines"
}

# lock_check PREFIX2 NAME STAGING LOGDIR
#   Rebuilds the environment from the two record files (conda explicit lock, then pip
#   --no-deps from the freeze list) and checks that `pip list --format=freeze` of the
#   rebuilt environment equals the recorded one.
lock_check() {
    local prefix2="$1" name="$2" staging="$3" logdir="$4"
    local conda="$EBROOTMINIFORGE3/bin/conda"
    "$conda" create -y -q -p "$prefix2" --file "$staging/$name.lock.txt" > "$logdir/$name.lockcheck.log" 2>&1
    "$prefix2/bin/python" -m pip install --no-deps -r "$staging/$name.pip-freeze.txt" >> "$logdir/$name.lockcheck.log" 2>&1
    "$prefix2/bin/python" -m pip list --format=freeze > "$TMPDIR/$name.lockcheck.freeze.txt"
    if diff "$staging/$name.pip-freeze.txt" "$TMPDIR/$name.lockcheck.freeze.txt" >> "$logdir/$name.lockcheck.log"; then
        echo "LOCKCHECK $name: rebuilt from lock + freeze, pip list identical" | tee -a "$logdir/$name.lockcheck.log"
    else
        echo "LOCKCHECK $name: pip list DIFFERS (see $name.lockcheck.log)" | tee -a "$logdir/$name.lockcheck.log"
    fi
}
