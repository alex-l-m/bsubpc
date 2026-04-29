{ config, pkgs, lib, ... }:

let
  # Plain Python interpreter only.
  # Project packages such as rdkit, numpy, pandas, etc. should be managed
  # per-project with uv, not installed globally through Nix.
  pythonForUv = pkgs.python312;

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
    glibc

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

  environment.sessionVariables = {
    # Let uv manage project environments, but make it use the nixpkgs Python
    # already installed in this VM rather than downloading generic Python
    # builds.
    UV_PYTHON_DOWNLOADS = "never";
    UV_PYTHON_PREFERENCE = "only-system";

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

      pkgs.codex
    ];

  system.stateVersion = "25.11";
}
