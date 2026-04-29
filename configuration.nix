{ config, pkgs, lib, ... }:

let
  # Plain Python interpreter only.
  # Project packages such as rdkit, numpy, pandas, etc. should be managed
  # per-project with uv, not installed globally through Nix.
  pythonForUv = pkgs.python312;

  # Runtime libraries that help generic Linux binaries / PyPI wheels run on
  # NixOS. This is plumbing, not Python dependency management.
  runtimeLibs = with pkgs; [
    stdenv.cc.cc.lib  # libstdc++.so.6, libgcc_s.so.1
    zlib
    bzip2
    xz
    openssl
    libffi
    sqlite
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

  # A very large amount of disk space to hold the CSD database
  virtualisation.diskSize = 64 * 1024; # MiB

  # Only like an eigth of my RAM; I don't think I'm doing anything
  # memory-constrained
  virtualisation.memorySize = 8 * 1024;    # MiB

  # Half my cores
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

  # Helps unpatched/generic Linux binaries and binary wheels find dynamic
  # libraries on NixOS.
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

    # Pragmatic convenience for a disposable VM. This helps some PyPI binary
    # wheels find libstdc++ and friends.
    LD_LIBRARY_PATH = lib.makeLibraryPath runtimeLibs;
  };

  environment.systemPackages =
    [
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
