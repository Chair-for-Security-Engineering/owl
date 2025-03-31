{ pkgs ? import <nixpkgs> {}}:

with pkgs; mkShell {
  NIX_ENFORCE_NO_NATIVE = 0;

  buildInputs = [
    gcc
    cmake
    (python3.withPackages(p: [ p.sympy ]))
  ];
}
