@echo off
REM pHinder - double-click launcher (Windows).
REM First run creates the 'pHinder' conda environment; later runs just open the app.
REM Requirement: install Miniforge once: https://conda-forge.org/download/
setlocal
set ENV_NAME=phinder
set HERE=%~dp0
set REPO=%HERE%..

set CONDA=
where conda >nul 2>nul && set CONDA=conda
if "%CONDA%"=="" (
  for %%C in ("%USERPROFILE%\miniforge3" "%USERPROFILE%\mambaforge" "%USERPROFILE%\miniconda3" "%USERPROFILE%\anaconda3") do (
    if exist "%%~C\Scripts\conda.exe" set CONDA="%%~C\Scripts\conda.exe"
  )
)
if "%CONDA%"=="" (
  echo Could not find conda on this PC.
  echo Please install Miniforge first: https://conda-forge.org/download/
  pause
  exit /b 1
)

%CONDA% env list | findstr /b /c:"%ENV_NAME% " >nul
if errorlevel 1 (
  echo First-time setup: creating the %ENV_NAME% environment...
  pushd "%REPO%"
  %CONDA% env create -f environment.yml
  set CREATE_ERR=%errorlevel%
  popd
  if not "%CREATE_ERR%"=="0" ( echo Setup did not finish. & pause & exit /b 1 )
)

REM --- keep this copy current -------------------------------------------------
REM Best-effort throughout: an offline PC, or a clone with local edits, still
REM launches on the code it already has.
set DOPULL=1
set NEEDS_ENV=
where git >nul 2>nul || set DOPULL=
if defined DOPULL git -C "%REPO%" rev-parse --is-inside-work-tree >nul 2>nul || set DOPULL=
if defined DOPULL git -C "%REPO%" remote get-url origin >nul 2>nul || set DOPULL=
if defined DOPULL for /f %%S in ('git -C "%REPO%" status --porcelain 2^>nul ^| find /c /v ""') do if not "%%S"=="0" set DOPULL=
REM A clone made with "--branch <tag>" sits on a detached HEAD. It can never
REM fast-forward, so it would stay on that release for ever, and silently.
REM Move it back onto the default branch when nothing can be lost.
if defined DOPULL git -C "%REPO%" symbolic-ref -q HEAD >nul 2>nul || call :unpin
if defined DOPULL echo Checking for updates...
if defined DOPULL for /f %%H in ('git -C "%REPO%" rev-parse HEAD 2^>nul') do set BEFORE=%%H
if defined DOPULL git -C "%REPO%" pull --ff-only --quiet 2>nul || echo   could not reach the server - launching the copy you have.
if defined DOPULL for /f %%H in ('git -C "%REPO%" rev-parse HEAD 2^>nul') do set AFTER=%%H
REM A new dependency is the only thing an editable install cannot pick up on its
REM own, so refresh the environment exactly when that pull touched the file.
if defined DOPULL if not "%BEFORE%"=="%AFTER%" for /f %%F in ('git -C "%REPO%" diff --name-only %BEFORE% %AFTER% 2^>nul ^| findstr /x "environment.yml"') do set NEEDS_ENV=1
if defined NEEDS_ENV echo Dependencies changed - updating the %ENV_NAME% environment...
if defined NEEDS_ENV pushd "%REPO%"
if defined NEEDS_ENV %CONDA% env update -f environment.yml
if defined NEEDS_ENV popd

echo Starting pHinder...
set PYTHONPATH=
cd /d "%USERPROFILE%"
%CONDA% run --no-capture-output -n %ENV_NAME% phinder-gui
exit /b

REM Move a tag-pinned clone back onto the default branch, but only when the
REM commit it sits on is already contained in that branch, so nothing is lost.
:unpin
set BRANCH=
for /f "tokens=3" %%B in ('git -C "%REPO%" remote show origin 2^>nul ^| findstr /c:"HEAD branch"') do set BRANCH=%%B
if "%BRANCH%"=="" set BRANCH=main
git -C "%REPO%" fetch --quiet origin %BRANCH% 2>nul || (set DOPULL=& goto :eof)
set HEADC=
for /f %%H in ('git -C "%REPO%" rev-parse HEAD 2^>nul') do set HEADC=%%H
git -C "%REPO%" merge-base --is-ancestor %HEADC% origin/%BRANCH% 2>nul || (
  echo This copy sits on a commit that is not part of '%BRANCH%' - leaving it alone.
  set DOPULL=
  goto :eof
)
echo This copy was pinned to a fixed release, which cannot receive updates.
echo Moving it onto '%BRANCH%' so it can...
git -C "%REPO%" checkout --quiet %BRANCH% 2>nul || set DOPULL=
goto :eof
