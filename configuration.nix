{ config, pkgs, lib, ... }:

let
  # Plain Python interpreter only.
  # Project packages such as rdkit, numpy, pandas, etc. should be managed
  # per-project with uv, not installed globally through Nix.
  pythonForUv = pkgs.python312;

  # The Grimme-lab xtb4stda + stda binaries are not in nixpkgs, but the
  # statically-linked Linux binaries on GitHub run as-is on NixOS. Wrap them
  # so XTB4STDAHOME points at the parameter files bundled in the source repo.
  xtb4stdaParams = pkgs.fetchFromGitHub {
    owner = "grimme-lab";
    repo = "xtb4stda";
    rev = "3b3d690d3726a1d16515a56cae14babf77bdcb5b";
    hash = "sha256-hDYgi6+QIJYmrsXDHWJdO2DCFrI91HoJ7oHsOSGLOY8=";
  };
  xtb4stdaBin = pkgs.fetchurl {
    url = "https://github.com/grimme-lab/xtb4stda/releases/download/v1.0/xtb4stda";
    hash = "sha256-fs5be4lZqKQGdwn23yNIDrMCxI5XogsddzZ8djNtMCs=";
  };
  stdaBin = pkgs.fetchurl {
    url = "https://github.com/grimme-lab/std2/releases/download/v1.6.3/stda_v1.6.3";
    hash = "sha256-lwHx0I553VYPRPuuYdEm6qzZROuPZ/ygggYhuQAzr4A=";
  };
  xtb4stda = pkgs.stdenvNoCC.mkDerivation {
    pname = "xtb4stda";
    version = "1.0-stda1.6.3";
    dontUnpack = true;
    nativeBuildInputs = [ pkgs.makeWrapper ];
    installPhase = ''
      mkdir -p $out/bin $out/share/xtb4stda
      install -m 0755 ${xtb4stdaBin} $out/bin/.xtb4stda-unwrapped
      install -m 0755 ${stdaBin} $out/bin/.stda-unwrapped
      install -m 0644 ${xtb4stdaParams}/.param_stda1.xtb $out/share/xtb4stda/
      install -m 0644 ${xtb4stdaParams}/.param_stda2.xtb $out/share/xtb4stda/
      install -m 0644 ${xtb4stdaParams}/.xtb4stdarc $out/share/xtb4stda/
      for gbsa in ${xtb4stdaParams}/.param_gbsa_*; do
        install -m 0644 "$gbsa" $out/share/xtb4stda/
      done
      # The xtb4stda binary uses a Fortran fixed-length string buffer for
      # XTB4STDAHOME (~80 chars), which truncates a /nix/store/.../share/xtb4stda
      # path. The wrappers default XTB4STDAHOME to /etc/xtb4stda — populated
      # below with environment.etc — so the runtime path stays short.
      makeWrapper $out/bin/.xtb4stda-unwrapped $out/bin/xtb4stda \
        --set-default XTB4STDAHOME /etc/xtb4stda
      makeWrapper $out/bin/.stda-unwrapped $out/bin/stda \
        --set-default XTB4STDAHOME /etc/xtb4stda
    '';
  };

  # CCDC install location used by the installer script.
  ccdcPrefix = "/opt/CCDC";
  ccdcPythonApiDir = "${ccdcPrefix}/ccdc-software/csd-python-api";

  # Runtime libraries for unpatched/generic Linux binaries through nix-ld.
  #
  # These are not Python dependency management. They are shared-library
  # plumbing for vendor binaries, bundled conda environments, and binary
  # wheels that expect a more FHS-like Linux system.
  runtimeLibs = with pkgs; [
    # C/C++ runtime.
    stdenv.cc.cc.lib  # libstdc++.so.6, libgcc_s.so.1
    gfortran.cc.lib   # libgfortran.so.5 for Fortran-backed extensions
    glibc

    # Linear algebra runtime for source-built tblite.
    blas
    lapack

    # Common compression / crypto / database / FFI libraries.
    zlib
    zstd
    bzip2
    xz
    openssl
    libffi
    sqlite
    expat

    # Needed by CCDC's bundled Qt networking stack:
    # provides libgssapi_krb5.so.2.
    krb5

    # Other common native-runtime libraries that often show up in
    # generic Linux binaries or binary Python wheels.
    curl
    ncurses
    util-linux  # libuuid, etc.

    # Qt / X11 / XCB runtime libraries needed by the CCDC installer and
    # potentially by CCDC's bundled Qt libraries.
    libxkbcommon
    libxcb
    libxcb-util
    libxcb-wm
    libxcb-image
    libxcb-keysyms
    libxcb-render-util
    xcb-util-cursor

    libx11
    libxext
    libxrender
    libxi
    libxfixes
    libxcursor
    libxrandr
    libxinerama
    libxau
    libxdmcp
    libsm
    libice

    # Font / GUI-adjacent libraries commonly needed by bundled Qt binaries.
    fontconfig
    freetype
    glib
    dbus
    libGL
  ];

  runtimeLibraryPath = lib.makeLibraryPath runtimeLibs;
  buildLibraryPath = lib.makeLibraryPath (with pkgs; [ blas lapack ]);
in
{
  # Same basic VM shape as the nix.dev tutorial.
  # https://github.com/NixOS/nix.dev/blob/master/source/tutorials/nixos/nixos-configuration-on-vm.md
  boot.loader.systemd-boot.enable = true;
  boot.loader.efi.canTouchEfiVariables = true;

  # Headless VM.
  virtualisation.graphics = false;

  # Temporarily needed for this bug:
  # https://github.com/NixOS/nixpkgs/issues/499166
  documentation.doc.enable = false;

  # CCDC's installer and generated wrappers use hardcoded shebangs such as
  # #!/bin/bash. envfs makes /bin/bash and similar paths resolve on NixOS
  # according to the process PATH.
  services.envfs.enable = true;

  # A very large amount of disk space to hold the CSD database.
  virtualisation.diskSize = 64 * 1024; # MiB

  # Only like an eighth of my RAM; I don't think I'm doing anything
  # memory-constrained.
  virtualisation.memorySize = 8 * 1024; # MiB

  # Half my cores.
  virtualisation.cores = 16;

  virtualisation.sharedDirectories.vmShare = {
    source = "/home/alexlm/vm-shared";
    target = "/mnt/host";
    securityModel = "none";
  };

  # This module provides QEMU-VM-specific options such as
  # virtualisation.forwardPorts.
  imports = [
    <nixpkgs/nixos/modules/virtualisation/qemu-vm.nix>
  ];

  # Run an SSH server inside the VM.
  services.openssh = {
    enable = true;

    settings = {
      PasswordAuthentication = false;
      KbdInteractiveAuthentication = false;
      PermitRootLogin = "no";
    };
  };

  # Make the VM's SSH port reachable from your host machine:
  #
  #   host localhost:2222  -->  VM port 22
  #
  # After the VM is running, connect with:
  #
  #   ssh -i ~/.ssh/id_alexlm -p 2222 alexlm@localhost
  virtualisation.forwardPorts = [
    {
      from = "host";
      host.port = 2222;
      guest.port = 22;
    }
  ];

  users.users.alexlm = {
    isNormalUser = true;
    description = "Alex";
    extraGroups = [ "wheel" ];

    openssh.authorizedKeys.keys = [
      "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIN/c+S1Osh+QIBAQBy3iWB3n134WrFcZhdzOROJuf29k alexlovesmolecules@gmail.com"
    ];
  };

  # Convenience for a disposable local VM: lets the wheel user run sudo
  # without needing a password that we would otherwise have to define.
  security.sudo.wheelNeedsPassword = false;

  # Helps unpatched/generic Linux binaries find dynamic libraries on NixOS.
  #
  # This is what lets the CCDC installer, activation binary, bundled conda
  # Python, and similar vendor binaries find libxkbcommon-x11.so.0,
  # libgssapi_krb5.so.2, libstdc++.so.6, etc.
  programs.nix-ld = {
    enable = true;
    libraries = runtimeLibs;
  };

  # xtb4stda needs its parameter files in a directory short enough to fit in
  # its Fortran path buffer. /etc/xtb4stda symlinks into the nix store copy.
  environment.etc."xtb4stda".source = "${xtb4stda}/share/xtb4stda";

  environment.sessionVariables = {
    # Let uv manage project environments, but make it use the nixpkgs Python
    # already installed in this VM rather than downloading generic Python
    # builds.
    UV_PYTHON_DOWNLOADS = "never";
    UV_PYTHON_PREFERENCE = "only-system";

    # programs.nix-ld owns NIX_LD_LIBRARY_PATH on current NixOS and points it at
    # /run/current-system/sw/share/nix-ld/lib. Binary Python wheels imported by
    # uv's venv still go through the normal dynamic loader path, so export the
    # same runtime closure as LD_LIBRARY_PATH for new shells.
    LD_LIBRARY_PATH = runtimeLibraryPath;

    # tblite is built from source by uv; Meson's cc.find_library needs these
    # libraries on the compiler search path.
    LIBRARY_PATH = buildLibraryPath;

    # CCDC Python API on a headless VM.
    #
    # This lets `import ccdc` work without trying to create a Qt QApplication.
    # Display-dependent API features still need a display or virtual display.
    CCDC_PYTHON_API_NO_QAPPLICATION = "1";
  };

  # Add CCDC's Python API wrapper directory to PATH.
  #
  # This gives you:
  #
  #   run_csd_python_api
  #
  # from any new login shell, without adding the bundled conda bin directory
  # directly to PATH.
  environment.extraInit = ''
    case ":$PATH:" in
      *:${ccdcPythonApiDir}:*) ;;
      *) export PATH="${ccdcPythonApiDir}:$PATH" ;;
    esac
  '';

  environment.systemPackages =
    [
      # Make sure `bash` is on PATH so envfs can resolve /bin/bash for
      # vendor scripts with hardcoded shebangs.
      pkgs.bashInteractive

      # Existing silliness.
      pkgs.cowsay
      pkgs.lolcat

      # Basic tools.
      pkgs.vim
      pkgs.git
      pkgs.curl
      pkgs.wget
      pkgs.unzip
      pkgs.file
      pkgs.ripgrep
      pkgs.fd

      # Binary/runtime diagnostics for future vendor-library issues.
      pkgs.binutils   # readelf, objdump
      pkgs.patchelf
      pkgs.strace

      # `magick montage` for compositing the per-molecule tiles emitted
      # by gap_grid.py into one PNG (see gap_grid.sh).
      pkgs.imagemagick

      # `xtb4stda` + `stda` binaries (Grimme-lab) used by stda_pipeline.sh.
      # See the let-binding above for how XTB4STDAHOME is wired up.
      xtb4stda

      # Build tools useful when Python/R packages need compilation.
      pkgs.pkg-config
      pkgs.gcc
      pkgs.gnumake
      pkgs.gfortran

      # Python tooling.
      #
      # uv manages project environments.
      # pythonForUv is only the base interpreter uv should use.
      pkgs.uv
      pythonForUv

      # R + the packages used by the analysis scripts (tidyverse for the bulk
      # of plotting, ggrepel for label placement in compare_stda.R, etc.).
      # `ggdark` is currently flagged broken in nixpkgs; compare_stda.R is
      # patched to drop the dark theme so we don't need it. Build the full
      # set with `nixos-rebuild test --max-jobs 1` to fit the small tmpfs
      # nix store on this VM.
      (pkgs.rWrapper.override {
        packages = with pkgs.rPackages; [
          tidyverse
          cowplot
          ggrepel
          robustbase
          tidymodels
          glue
          matsindf
        ];
      })

      pkgs.codex
    ];

  system.stateVersion = "25.11";
}
