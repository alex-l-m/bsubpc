uv sync

nix \
    --extra-experimental-features nix-command \
    --extra-experimental-features flakes \
    shell nixpkgs#patchelf \
    --command patchelf \
    --add-needed \
    /nix/store/fsri8ivf3jgynw975x8k57b5arznkav1-gfortran-15.2.0-lib/lib/libgfortran.so \
    .venv/lib/python3.12/site-packages/tblite/_libtblite.cpython-312-x86_64-linux-gnu.so
