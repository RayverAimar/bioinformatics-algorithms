Write-Host "[SCRIPT] Compiling global_alignment.cpp..."
& g++ -o global_alignment.exe global_alignment.cpp

if ($LASTEXITCODE -ne 0) {
    Write-Host "Compilation failed."
    pause
    exit 1
}

Write-Host "`n[SCRIPT] Compilation successful. Running global_alignment executable..."
Start-Process -FilePath ".\global_alignment.exe" -NoNewWindow -Wait

Write-Host "`n[SCRIPT] Running main.py..."
python.exe dot_plot_maker.py

pause