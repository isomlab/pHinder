#!/bin/bash
# pHinder — double-click launcher (macOS).
#
# First run: creates the 'pHinder' conda environment from environment.yml
# (Python + the app + its dependencies); this can take a few minutes.
# Every run after that: just opens the app.
#
# Requirement: install Miniforge once (a normal clickable installer):
#   https://conda-forge.org/download/

ENV_NAME="phinder"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"   # this script lives in <repo>/launchers/

pause_and_exit() {
    echo
    read -r -p "Press Return to close this window…" _
    exit "${1:-1}"
}

find_conda() {
    local c
    for c in "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" \
             "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
             "/opt/homebrew/Caskroom/miniforge/base/bin/conda" \
             "$(command -v conda 2>/dev/null)"; do
        if [ -n "$c" ] && [ -x "$c" ]; then
            echo "$c"; return 0
        fi
    done
    return 1
}

CONDA="$(find_conda)" || {
    echo "Could not find conda on this Mac."
    echo "Please install Miniforge first (clickable installer):"
    echo "    https://conda-forge.org/download/"
    pause_and_exit 1
}

# Create the environment the first time only.
if ! "$CONDA" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "First-time setup: creating the '$ENV_NAME' environment (a few minutes)…"
    echo
    ( cd "$REPO" && "$CONDA" env create -f environment.yml ) || {
        echo
        echo "Setup did not finish. Please see the messages above."
        pause_and_exit 1
    }
    echo
    echo "Setup complete."
fi

# --- keep this copy current -------------------------------------------------
# Best-effort throughout: an offline laptop, or a clone with local edits, still
# launches on the code it already has.

update_repo() {
    command -v git >/dev/null 2>&1 || return 0
    git -C "$REPO" rev-parse --is-inside-work-tree >/dev/null 2>&1 || return 0
    git -C "$REPO" remote get-url origin >/dev/null 2>&1 || return 0
    if [ -n "$(git -C "$REPO" status --porcelain 2>/dev/null)" ]; then
        echo "This copy has local changes — skipping update."
        return 0
    fi
    # A clone made with "--branch <tag>" sits on a detached HEAD. It can never
    # fast-forward, so it would stay on that release for ever — and silently,
    # which is worse. Move it back onto the default branch, but only when
    # nothing can be lost: the tree is clean (checked above) and the commit it
    # is sitting on is already contained in that branch.
    if ! git -C "$REPO" symbolic-ref -q HEAD >/dev/null 2>&1; then
        local branch head
        branch="$(git -C "$REPO" remote show origin 2>/dev/null | awk '/HEAD branch/ {print $NF}')"
        [ -n "$branch" ] || branch="main"
        if ! git -C "$REPO" fetch --quiet origin "$branch" 2>/dev/null; then
            echo "  could not reach the server — launching the copy you have."
            return 0
        fi
        head="$(git -C "$REPO" rev-parse HEAD 2>/dev/null)"
        if ! git -C "$REPO" merge-base --is-ancestor "$head" "origin/$branch" 2>/dev/null; then
            echo "This copy sits on a commit that is not part of '$branch' — leaving it alone."
            return 0
        fi
        echo "This copy was pinned to a fixed release, which cannot receive updates."
        echo "Moving it onto '$branch' so it can…"
        if ! git -C "$REPO" checkout --quiet "$branch" 2>/dev/null; then
            echo "  could not switch — launching the copy you have."
            return 0
        fi
    fi
    echo "Checking for updates…"
    local before after
    before="$(git -C "$REPO" rev-parse HEAD 2>/dev/null)"
    if ! git -C "$REPO" pull --ff-only --quiet 2>/dev/null; then
        echo "  could not reach the server — launching the copy you have."
        return 0
    fi
    after="$(git -C "$REPO" rev-parse HEAD 2>/dev/null)"
    if [ "$before" = "$after" ]; then echo "  already up to date."; else echo "  updated."; fi
}

# The app is installed editable, so new code needs no reinstall — a new
# dependency does. Rebuild only when environment.yml is newer than the env,
# which also covers a pull done by hand outside this launcher.
update_env() {
    local prefix yml
    yml="$REPO/environment.yml"
    [ -f "$yml" ] || return 0
    prefix="$("$CONDA" env list | awk -v n="$ENV_NAME" '$1 == n {print $NF}')"
    [ -n "$prefix" ] && [ -f "$prefix/conda-meta/history" ] || return 0
    if [ "$yml" -nt "$prefix/conda-meta/history" ]; then
        echo "Dependencies changed — updating the '$ENV_NAME' environment…"
        ( cd "$REPO" && "$CONDA" env update -f environment.yml ) \
            || echo "  update failed — launching on the environment you have."
    fi
}

update_repo
update_env

echo "Starting pHinder…"
# Isolate from the user's Python environment: PYTHONPATH is cleared (entries there
# take precedence over the env's site-packages and can shadow the app's package),
# and we run from a neutral directory (conda run also puts the current dir on
# sys.path). The created environment is self-contained.
cd "$HOME" || cd /
exec env -u PYTHONPATH "$CONDA" run --no-capture-output -n "$ENV_NAME" phinder-gui
