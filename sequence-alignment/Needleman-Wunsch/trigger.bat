@echo off
echo [SCRIPT] Compiling global_alignment.cpp...
g++ -o global_alignment.exe global_alignment.cpp

if %errorlevel% neq 0 (
  echo Compilation failed.
  pause
  exit /b 1
)

echo.
echo [SCRIPT] Compilation successful. Running global_alignment executable...
global_alignment.exe

echo.
echo [SCRIPT] Running dot_plot_maker.py...
python.exe dot_plot_maker.py

pause